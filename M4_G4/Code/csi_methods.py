# M4 NOTE: functions below the marker "# ---- M4 additions ----" are new in M4; everything above it is the unmodified M3 code
# (see Code/original_provided/csi_methods.py). Nothing above the marker was changed.
"""Phase-correction methods of Ratnam et al. (IEEE TWC 2024) written for REAL CSI whose subcarrier indices may have gaps
(e.g. the DC/null tones of 802.11ax RU996). Self-contained (numpy only). Index vector `kidx` gives the integer subcarrier index of
each column; the subcarrier frequency offset is f_k = kidx / Ts.

Sign convention (derived from the paper's Eq. 20 and Eq. 3):  h_bar = h * exp(-j 2 pi f tau) * exp(-j psi);
the corrected CSI is  h_hat = h_bar * exp(+j 2 pi f tau_hat) * exp(+j psi_hat).
Equation numbers refer to the paper (arXiv 2307.12126 v2, accepted version)."""
import numpy as np

TS_HE = 12.8e-6                                        # 802.11ax HE-LTF 4x symbol without CP: 78.125 kHz spacing
HE80_KIDX = np.r_[np.arange(-500, -2), np.arange(3, 501)]   # ASSUMED tone map of RU996 (996 tones). CONFIRM against the dataset README / SHARPax code.

def blocks(kidx):
    """List of index arrays of contiguous runs of kidx (split where the index jumps by more than 1)."""
    cut = np.where(np.diff(kidx) != 1)[0] + 1
    return np.split(np.arange(len(kidx)), cut)

def gain_eq6(X):                                       # Eq.(6): divide each frame by sqrt of its mean power
    return X / np.sqrt(np.mean(np.abs(X) ** 2, axis=1, keepdims=True))

def apply_phase(X, tau, psi, kidx, Ts=TS_HE):          # Eq.(3)
    f = kidx / Ts
    return X * np.exp(2j * np.pi * f[None, :] * tau[:, None]) * np.exp(1j * psi)[:, None]

def eq22(X, kidx, Ts=TS_HE):
    """Eq.(22a),(22b). Only pairs of ADJACENT indices (difference 1) are used in the sum, so the DC gap is skipped."""
    j = np.where(np.diff(kidx) == 1)[0]
    s = np.sum(X[:, j] * np.conj(X[:, j + 1]), axis=1)
    tau = Ts / (2 * np.pi) * np.angle(s)
    f = kidx / Ts
    psi = -np.angle(np.sum(X * np.exp(2j * np.pi * f[None, :] * tau[:, None]), axis=1))
    return tau, psi

def eq21(X, kidx, Ts=TS_HE, noff=range(-3, 4)):
    """Eq.(21): least-squares line fit to the unwrapped phase, per frame. Unwrapping is done inside each contiguous block; for a second
    block the unknown 2*pi offset n is searched (-3..3) and the offset with the smallest residual is kept. Phase model: -2 pi f tau - psi."""
    f = kidx / Ts; ph = np.angle(X); P = X.shape[0]
    bl = blocks(kidx)
    U = np.empty_like(ph)
    for b in bl: U[:, b] = np.unwrap(ph[:, b], axis=1)
    A = np.stack([-2 * np.pi * f, -np.ones_like(f)], axis=1)        # params [tau, psi]
    Ap = np.linalg.pinv(A)
    best, bestsse = None, None
    opts = noff if len(bl) == 2 else [0]
    for n in opts:
        Y = U.copy()
        if len(bl) == 2: Y[:, bl[1]] += 2 * np.pi * n
        th = Y @ Ap.T
        sse = np.sum((Y - th @ A.T) ** 2, axis=1)
        if best is None: best, bestsse = th, sse
        else:
            m = sse < bestsse; best = np.where(m[:, None], th, best); bestsse = np.where(m, sse, bestsse)
    return best[:, 0], best[:, 1]

def algo4_eq25(X, kidx, Ts=TS_HE, b_ref=None, thr=0.1):
    """Algorithm 4 with Eq.(25),(26). X must be gain-corrected (mean power 1 per frame so that |b|^2 is comparable to the 0.1 threshold).
    b_ref: optional static reference (complex, length K) corrected-by-construction; if None it is estimated from X via Eq.(23)
    (coarse tau, psi from Eq.22). Returns tau_hat, psi_hat."""
    f = kidx / Ts
    tc, pc = eq22(X, kidx, Ts)
    if b_ref is None:
        b_ref = apply_phase(X, tc, pc, kidx, Ts).mean(axis=0)                                  # Eq.(23)
    kb = np.where(np.abs(b_ref) ** 2 > thr)[0]                                                 # K_bar
    om = np.conj(X[:, kb]) * b_ref[kb][None, :] * np.exp(-2j * np.pi * f[kb][None, :] * tc[:, None])
    # robust unwrapping, Eq.(26): local sum over the 3 preceding and 3 succeeding indices of K_bar
    S = np.apply_along_axis(lambda r: np.convolve(r, np.ones(7), mode="same"), 1, om)
    u = np.unwrap(np.angle(S), axis=1)
    phi = np.mod(np.angle(om) - u + np.pi, 2 * np.pi) - np.pi + u
    w = np.abs(om); fp = 2 * np.pi * f[kb][None, :]
    Sw, Sf, Sff = w.sum(1), (w * fp).sum(1), (w * fp * fp).sum(1)
    Sp, Sfp = (w * phi).sum(1), (w * fp * phi).sum(1)
    den = Sw * Sff - Sf ** 2
    dtau = (Sw * Sfp - Sf * Sp) / den                                                          # weighted LS of Eq.(25)
    psi = (Sff * Sp - Sf * Sfp) / den
    return tc + dtau, psi

def correct(method, X, kidx, Ts=TS_HE, b_ref=None):
    """X: raw complex CSI (P,K). Returns CSI after Eq.(6) gain normalisation and the chosen phase method.
    method in {'raw_nogain','A0','A1_eq21','A2_eq22','A3_algo4'}"""
    if method == "raw_nogain": return X
    G = gain_eq6(X)
    if method == "A0": return G
    if method == "A1_eq21": t, p = eq21(G, kidx, Ts)
    elif method == "A2_eq22": t, p = eq22(G, kidx, Ts)
    elif method == "A3_algo4": t, p = algo4_eq25(G, kidx, Ts, b_ref=b_ref)
    else: raise ValueError(method)
    return apply_phase(G, t, p, kidx, Ts)

def dynamic_power_fraction_db(X):
    """KPI (mentor's metric): share of CSI power that varies over time around its own mean, in dB (lower is better).
    V = sum_p sum_k |x[p,k] - mean_p x[:,k]|^2 / sum_p sum_k |x[p,k]|^2 (complex mean, so phase wobble counts as variation)."""
    mu = X.mean(axis=0, keepdims=True)
    return 10 * np.log10(np.sum(np.abs(X - mu) ** 2) / np.sum(np.abs(X) ** 2))

def phase_circ_std_deg(X, thr=0.1):
    """Secondary KPI: circular std (degrees) of the phase of each subcarrier around its temporal mean phase, median over subcarriers
    whose mean |.|^2 exceeds thr (after Eq.6 normalisation). Lower is better."""
    mu = X.mean(axis=0)
    m = np.abs(mu) ** 2 > thr * np.mean(np.abs(X) ** 2)
    R = np.abs(np.mean(np.exp(1j * (np.angle(X[:, m]) - np.angle(mu[m])[None, :])), axis=0))
    return float(np.degrees(np.median(np.sqrt(-2 * np.log(np.clip(R, 1e-9, 1))))))

# ---- M4 additions ----
# Tone map of the Brescia 20 MHz dataset (HE RU242, 256-point FFT, fft-shifted, 6 low guard + 3 DC + 5 high guard tones empty):
# MAT row r (0..255) is subcarrier index k = r - 128; usable rows 6..126 and 130..250  ->  k = -122..-2 and 2..122 (242 tones).
HE20_ROWS = np.r_[np.arange(6, 127), np.arange(130, 251)]
HE20_KIDX = HE20_ROWS - 128

def complex_consistency(X):
    """Complex-CSI consistency (phase-sensitive, higher is better, 1 = identical up to a real gain, ~0 = random common phase).
    C_mean: mean over frames of Re<x_p, mu>/(|x_p||mu|), mu = temporal mean.  C_lag1: mean over consecutive frames of
    Re<x_p, x_{p+1}>/(|x_p||x_{p+1}|)."""
    mu = X.mean(axis=0)
    nx = np.linalg.norm(X, axis=1)
    c_mean = np.mean(np.real(X @ np.conj(mu)) / (nx * np.linalg.norm(mu) + 1e-30))
    c_lag1 = np.mean(np.real(np.sum(X[:-1] * np.conj(X[1:]), axis=1)) / (nx[:-1] * nx[1:] + 1e-30))
    return float(c_mean), float(c_lag1)

def eq21_diagnostics(X_gain, kidx, Ts=TS_HE):
    """Per-frame Eq.(21) fit residual (rms of unwrapped-phase line fit, rad) - used only as a diagnostic of unwrapping quality."""
    t, p = eq21(X_gain, kidx, Ts)
    f = kidx / Ts
    ph = np.angle(apply_phase(X_gain, t, p, kidx, Ts))
    return t, p, np.sqrt(np.mean(ph ** 2, axis=1))

def static_tone_mask(Y_ref, thr=0.1):
    """Tones whose temporal-mean power (of a reference-corrected window, e.g. Algorithm 4 output) exceeds thr x mean power.
    Used so that every method's phase spread is computed on the SAME tones (for raw/A0 the temporal mean collapses, so a per-method mask is empty)."""
    return np.abs(Y_ref.mean(axis=0)) ** 2 > thr * np.mean(np.abs(Y_ref) ** 2)

def phase_circ_std_deg_masked(X, mask):
    """Circular std (deg) of subcarrier phase around its temporal mean phase, median over the given tone mask. If the temporal mean of a
    tone is ~0 its mean phase is meaningless: that is exactly what the KPI should penalise, so the value is computed, not skipped."""
    mu = X.mean(axis=0)
    R = np.abs(np.mean(np.exp(1j * (np.angle(X[:, mask]) - np.angle(mu[mask])[None, :])), axis=0))
    return float(np.degrees(np.median(np.sqrt(-2 * np.log(np.clip(R, 1e-9, 1))))))

def eq21_wide(X, kidx, Ts=TS_HE):
    """Eq.(21) exactly as in the M3 code but with the unknown 2*pi offset of the second tone block searched over -8..8 instead of -3..3
    (M4 finding: frames with a large timing error need |n| > 3; the M3 range clips them)."""
    return eq21(X, kidx, Ts, noff=range(-8, 9))

def correct_m4(method, X, kidx, Ts=TS_HE, b_ref=None):
    """Like correct(), plus method 'A1w_eq21_wide' (Eq.6 + Eq.21 with the widened offset search)."""
    if method == "A1w_eq21_wide":
        G = gain_eq6(X); t, p = eq21_wide(G, kidx, Ts); return apply_phase(G, t, p, kidx, Ts)
    return correct(method, X, kidx, Ts, b_ref)
