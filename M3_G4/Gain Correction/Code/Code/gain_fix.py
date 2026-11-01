"""
gain_fix.py -- M3 (Jenil): investigation and correction of Algorithm 3 of
V. V. Ratnam et al., "Optimal preprocessing of WiFi CSI for sensing applications", IEEE TWC 2024
(Sec. III-D, Eqs. 14-19, Algorithms 2 and 3).

csi_sim.py (M2) is NOT modified; it is imported read-only.

Public API
----------
gain_algo3(ht, variant="fixed"|"original", fixes=None, ...)      -> (gain[P], info)
gain_algo3_from_power(G, ...)                                     -> (gain or None, info)
gain_eq8(ht)                     baseline (8), needed for the Figure 6 comparison
chi_metric_fast(hh, sim)         Eq. (29); numerically identical to csi_sim.chi_metric (cached tau grid)
gain_rms_error_db(g_hat, g)      direct gain-tracking error (dB, constant offset removed)

variant="original" reproduces csi_sim.gain_algo3 exactly (checked in tests/test_gain.py) and additionally
returns the selected lambda and the full per-lambda table.

variant="fixed" (default) applies FIX-1 ... FIX-5:
  FIX-1 "grid"     lambda grid = multiples of 0.01 dB inside the paper's range [0.05, 1] x lambda_max.
                   (paper / M2: 20 points spaced 0.05*lambda_max = 0.14-0.21 dB; 0.5 dB is never on it.)
  FIX-2 "rayleigh" acceptance test = Rayleigh test of circular non-uniformity, P*m^2 > ln(1/alpha),
                   alpha = 1e-3, m = |mean_p exp(j 2 pi Gh_p / lambda)| (the statistic of Eq. 18).
                   (paper / M2: mean(Gh^2) <= lambda^2/24, which rejects the true lambda as soon as the
                   residual std exceeds ~0.2 lambda = 0.1 dB, although Eq. 18 is still accurate there.)
  FIX-3 "tiebreak" objective ties (within 1e-9 dB^2) -> largest lambda (lambda and lambda/n correct the
                   AGC identically without noise; Remark 1 favours larger lambda).
  FIX-4 "round"    Eq. (16) rounding floor(x + 1/2) (round-half-up) instead of np.round (half-to-even).
                   The "-1/2" outside the floor in the printed Eq. (16) is a typo (see tests).
  FIX-5            explicit fallback flag (info["fallback"]) when no lambda is accepted (M2: silent Eq. 6).
Optional, NOT in the default fix set (evaluated as ablations, see docs/gain_equation_mapping.md):
  "bonf"  Bonferroni correction of FIX-2 over the grid (alpha / number of lambdas).
  "loo"   leave-one-out moving average inside Algorithm 2.
Algorithm 2 itself (LPF of Xi, one-sided width 6/Trep, unwrapped angle) is unchanged in all variants.
"""
import numpy as np
from scipy.stats import norm
from sklearn.cluster import DBSCAN
import csi_sim as c

TREP = c.TREP
ALL_FIXES = frozenset({"grid", "rayleigh", "tiebreak", "round"})
_ZZ = np.arange(-30, 31)                       # same truncation |z| <= 30 as csi_sim._D


# ------------------------------------------------------------------ building blocks
def lpf_ma_fast(x, half_width, axis=-1):
    """Vectorised csi_sim.lpf_ma (moving average, window truncated at the edges). Identical output."""
    x = np.moveaxis(np.asarray(x), axis, -1); n = x.shape[-1]
    cs = np.concatenate([np.zeros(x.shape[:-1] + (1,), x.dtype), np.cumsum(x, -1)], -1)
    i = np.arange(n); lo = np.maximum(0, i - half_width); hi = np.minimum(n, i + half_width + 1)
    return np.moveaxis((cs[..., hi] - cs[..., lo]) / (hi - lo), -1, axis)


def lpf_ma_loo(x, half_width, axis=-1):
    """Leave-one-out version (sample i excluded from its own window). Optional ablation 'loo'."""
    x = np.moveaxis(np.asarray(x), axis, -1); n = x.shape[-1]
    cs = np.concatenate([np.zeros(x.shape[:-1] + (1,), x.dtype), np.cumsum(x, -1)], -1)
    i = np.arange(n); lo = np.maximum(0, i - half_width); hi = np.minimum(n, i + half_width + 1)
    return np.moveaxis((cs[..., hi] - cs[..., lo] - x) / (hi - lo - 1), -1, axis)


def D(x):
    """Eq. (19): D(x) = sum_z [Q((z-1/2)x) - Q((z+1/2)x)] z^2 ; vectorised over x."""
    x = np.atleast_1d(np.asarray(x, float))[:, None]
    return np.sum((norm.sf((_ZZ - 0.5) * x) - norm.sf((_ZZ + 0.5) * x)) * _ZZ ** 2, axis=1)


def csi_power_db(ht):
    """Eq. (7): observed CSI power per frame in dB."""
    return 10 * np.log10(np.mean(np.abs(ht) ** 2, axis=1))


def algo2_g1(G, lam, hw, loo=False):
    """Algorithm 2 for one or many lambdas (lam scalar or array -> output shape (..., P))."""
    lam = np.asarray(lam, float)[..., None]
    Xi = np.exp(2j * np.pi * np.asarray(G)[None, :] / lam) if lam.ndim > 1 else np.exp(2j * np.pi * G / lam)  # Eq. (14)
    Xb = lpf_ma_loo(Xi, hw) if loo else lpf_ma_fast(Xi.real, hw) + 1j * lpf_ma_fast(Xi.imag, hw)
    return np.unwrap(np.angle(Xb), axis=-1) * lam / (2 * np.pi)


# ------------------------------------------------------------------ Algorithm 3
def gain_algo3_from_power(G, variant="fixed", fixes=None, hw=None, step_db=0.01, alpha=1e-3,
                          eq16_literal=False):
    """Algorithm 3 on the observed power series G (dB). Returns (gain or None, info)."""
    hw = hw or int(round(6 / TREP))                           # moving average, one-sided width 6/Trep
    fixes = frozenset((set() if variant == "original" else ALL_FIXES) if fixes is None else fixes)
    G = np.asarray(G, float); P = len(G)
    lam_max = 1.5 * (G.max() - G.min()) + 1e-9               # Algorithm 3, line 2
    if "grid" in fixes:                                      # FIX-1
        n0 = max(int(np.ceil(0.05 * lam_max / step_db)), 1); n1 = int(np.floor(lam_max / step_db))
        grid = np.round(np.arange(n0, n1 + 1) * step_db, 10)
        if grid.size == 0: grid = np.array([lam_max])
    else:
        grid = np.arange(0.05, 1.0001, 0.05) * lam_max         # (0.05:0.05:1) x lambda_max
    g1 = algo2_g1(G, grid, hw, loo="loo" in fixes)           # Algorithm 2, all lambdas at once
    x = (G[None, :] - g1) / grid[:, None]
    z = np.floor(x + 0.5) if "round" in fixes else np.round(x)  # Eq. (16) (nearest integer)
    if eq16_literal: z = z - 0.5                              # printed Eq. (16) taken literally (typo test)
    g2 = grid[:, None] * z
    Gh = G[None, :] - g1 - g2                                 # Eq. (17)
    v = np.mean(Gh ** 2, axis=1)
    m = np.abs(np.mean(np.exp(2j * np.pi * Gh / grid[:, None]), axis=1))
    if "rayleigh" in fixes:                                   # FIX-2 (+ optional Bonferroni)
        a = alpha / len(grid) if "bonf" in fixes else alpha
        acc = P * m ** 2 > np.log(1 / a)
    else:
        acc = v <= grid ** 2 / 24                             # Algorithm 3 test as printed
    s2 = -(grid ** 2) / (2 * np.pi ** 2) * np.log(np.maximum(m, 1e-12))      # Eq. (18)
    obj = np.full(len(grid), np.inf)
    if acc.any():
        obj[acc] = s2[acc] + grid[acc] ** 2 * D(grid[acc] / np.sqrt(np.maximum(s2[acc], 1e-12)))  # Obj, Eq. (19)
    if np.isfinite(obj).any():
        if "tiebreak" in fixes:                               # FIX-3
            best = np.where(obj <= obj.min() + 1e-9)[0].max()
        else:
            best = int(np.argmin(obj))                        # first minimum, as the M2 loop
    else:
        best = None
    table = np.column_stack([grid, v, grid ** 2 / 24, m, acc, np.where(acc, s2, np.nan), obj])
    info = dict(fixes=sorted(fixes), lam_max=lam_max, grid=grid, table=table,
                fallback=best is None, lam_hat=np.nan if best is None else float(grid[best]))
    if best is None:
        return None, info
    info["g1"], info["g2"] = g1[best], g2[best]
    return 10 ** ((g1[best] + g2[best]) / 20), info


def gain_algo3(ht, variant="fixed", **kw):
    gain, info = gain_algo3_from_power(csi_power_db(ht), variant=variant, **kw)
    if gain is None:                                          # FIX-5 (same Eq. 6 fallback as M2, but flagged)
        gain = c.gain_eq6(ht)
    return gain, info


# ------------------------------------------------------------------ Eq. (8) baseline (Figure 6 comparison)
def gain_eq8(ht, eps=0.15):
    """Baseline (8): DBSCAN on the observed power levels (eps = 0.15, min-points = 1), cluster mean."""
    G = csi_power_db(ht)
    lab = DBSCAN(eps=eps, min_samples=1).fit_predict(G.reshape(-1, 1))
    out = np.empty_like(G)
    for l in np.unique(lab):
        out[lab == l] = G[lab == l].mean()
    return 10 ** (out / 20)


# ------------------------------------------------------------------ metrics
_TAUS = np.linspace(-4e-7, 4e-7, 1601); _EXP = {}
def chi_metric_fast(hh, sim):
    """Eq. (29) exactly as csi_sim.chi_metric, with the 1601 x K steering matrix cached per K."""
    b, d, f, K = sim["b"], sim["d"], sim["f"], hh.shape[1]
    P = hh.shape[0]; gamma = sim["gamma"]
    if K not in _EXP: _EXP[K] = np.exp(2j * np.pi * f[None, :] * _TAUS[:, None])
    bh = hh.mean(0)
    tau = _TAUS[np.argmax(np.abs(_EXP[K] @ (b * np.conj(bh))))]
    num = np.abs(np.sum(np.conj(hh - bh[None, :]) * d * np.exp(2j * np.pi * f[None, :] * tau))) ** 2
    den = (1 - gamma) * K * P * np.sum(np.abs(hh - bh[None, :]) ** 2)
    return float(num / den)


def gain_rms_error_db(g_hat, g_true):
    """RMS gain-tracking error in dB after removing the constant offset (a common scale is irrelevant)."""
    e = 20 * np.log10(g_hat) - 20 * np.log10(g_true); e = e - e.mean()
    return float(np.sqrt(np.mean(e ** 2)))


def sigma_gamma_db(K, gamma, corrected=True):
    """Std (dB) of the true CSI power Gamma_p. Paper Lemma 1: 100(1-g^2)/K. Corrected (ln10): 18.86(1-g^2)/K."""
    k = (10 / np.log(10)) ** 2 if corrected else 100.0
    return float(np.sqrt(k * (1 - gamma ** 2) / K))
