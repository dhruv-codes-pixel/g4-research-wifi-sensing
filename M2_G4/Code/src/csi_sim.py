"""
M2 artifact: simulation + reproduction of Ratnam et al., "Optimal preprocessing of WiFi CSI
for sensing applications", IEEE TWC 2024 (arXiv:2307.12126), Sec. V-B style experiments.

CSI model (paper eq. 2):  h~[p,k] = g_p * h[p,k] * exp(-j 2 pi f_k tau_p) * exp(-j psi_p)
Implemented here:
  gain correction : baseline Eq.(6), Algorithm 1 (DBSCAN on power increments + 0.1 Hz LPF)
  phase correction: baseline Eq.(21), baseline Eq.(22), Algorithm 4 + Eq.(25)
Metric: correlation chi (Eq.29) and SNR = chi/(1-chi)   (chi as defined in the paper is a squared correlation)
"""
import numpy as np
from sklearn.cluster import DBSCAN

TS = 3.2e-6          # OFDM symbol duration (s) -> subcarrier spacing 312.5 kHz
TREP = 0.1           # inter-frame interval (s)

# ---------------------------------------------------------------- simulator
def make_static(K, gamma, rng, max_delay=200e-9, n_taps=6):
    """Static channel b_k: random exponential-PDP multipath (SIMPLIFICATION of 802.11ax Model-C)."""
    f = np.arange(K) / TS
    delays = np.sort(rng.uniform(0, max_delay, n_taps)); delays[0] = 0.0
    amps = (rng.standard_normal(n_taps) + 1j * rng.standard_normal(n_taps)) * np.exp(-delays / 80e-9)
    amps[0] *= 3.0                                   # stronger first path
    b = (amps[None, :] * np.exp(-2j * np.pi * f[:, None] * delays[None, :])).sum(1)
    theta = np.angle(np.sum(b[:-1] * np.conj(b[1:])))  # enforce  angle(sum b_k b*_{k+1}) = 0
    b = b * np.exp(1j * theta * np.arange(K))
    b *= np.sqrt(gamma * K / np.sum(np.abs(b) ** 2))   # sum|b|^2 / K = gamma
    return b

def simulate(K=64, P=300, gamma=0.9, tau_max=1e-7, rng=None, agc=(-0.5, 0, 0.5), agc_p=(0.2, 0.6, 0.2)):
    rng = rng or np.random.default_rng()
    f = np.arange(K) / TS
    b = make_static(K, gamma, rng)
    d = np.sqrt((1 - gamma) / 2) * (rng.standard_normal((P, K)) + 1j * rng.standard_normal((P, K)))
    h = b[None, :] + d
    # large-scale gain: Gaussian process, std 0.2 dB, band-limited to <0.1 Hz (moving-average lowpass)
    w = rng.standard_normal(P + 120); ker = np.ones(21) / 21
    g1 = np.convolve(w, ker, mode="same")[60:60 + P]; g1 = 0.2 * g1 / g1.std()
    g2 = rng.choice(agc, size=P, p=agc_p)             # AGC gain (dB), discrete
    g = 10 ** ((g1 + g2) / 20)
    tau = rng.uniform(0, tau_max, P)
    psi = rng.uniform(-np.pi, np.pi, P)
    ht = g[:, None] * h * np.exp(-2j * np.pi * f[None, :] * tau[:, None]) * np.exp(-1j * psi)[:, None]
    return dict(ht=ht, h=h, b=b, d=d, g=g, tau=tau, psi=psi, f=f, gamma=gamma)

# ---------------------------------------------------------------- gain correction
def gain_eq6(ht):
    return np.sqrt(np.mean(np.abs(ht) ** 2, axis=1))

def lpf_ma(x, half_width):
    n = len(x); out = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - half_width), min(n, i + half_width + 1)
        out[i] = x[lo:hi].mean()
    return out

def gain_algo1(ht, eps=0.2):
    P = ht.shape[0]
    G = 10 * np.log10(np.mean(np.abs(ht) ** 2, axis=1))
    dG = np.zeros(P); dG[1:] = G[1:] - G[:-1]
    lab = DBSCAN(eps=eps, min_samples=1).fit_predict(dG[1:].reshape(-1, 1))
    cm = {l: dG[1:][lab == l].mean() for l in np.unique(lab)}
    g2 = np.zeros(P)
    for p in range(1, P):
        g2[p] = g2[p - 1] + cm[lab[p - 1]]
    g1 = lpf_ma(G - g2, half_width=int(round(6 / TREP)))     # 0.1 Hz moving-average LPF
    return 10 ** ((g1 + g2) / 20)


from scipy.stats import norm
def _D(x, zmax=30):
    z = np.arange(-zmax, zmax + 1)
    return np.sum((norm.sf((z - 0.5) * x) - norm.sf((z + 0.5) * x)) * z ** 2)

def gain_algo3(ht, hw=None):
    """Algorithm 3: ML gain estimation for uniformly spaced AGC grid (unknown step lambda, line search)."""
    hw = hw or int(round(6 / TREP))
    G = 10 * np.log10(np.mean(np.abs(ht) ** 2, axis=1))
    lam_max = 1.5 * (G.max() - G.min()) + 1e-9
    best = (np.inf, None)
    for lam in np.arange(0.05, 1.0001, 0.05) * lam_max:
        Xi = np.exp(2j * np.pi * G / lam)
        Xb = lpf_ma(Xi.real, hw) + 1j * lpf_ma(Xi.imag, hw)
        g1 = np.unwrap(np.angle(Xb)) * lam / (2 * np.pi)
        g2 = lam * np.round((G - g1) / lam)
        Gh = G - g1 - g2
        if np.mean(Gh ** 2) > lam ** 2 / 24:
            continue
        m = np.abs(np.mean(np.exp(2j * np.pi * Gh / lam)))
        sig2 = -(lam ** 2) / (2 * np.pi ** 2) * np.log(max(m, 1e-12))
        obj = sig2 + lam ** 2 * _D(lam / np.sqrt(max(sig2, 1e-12)))
        if obj < best[0]:
            best = (obj, 10 ** ((g1 + g2) / 20))
    return best[1] if best[1] is not None else gain_eq6(ht)

# ---------------------------------------------------------------- phase correction
def apply_phase(hb, f, tau_hat, psi_hat):                      # paper Eq.(3) without gain term
    return hb * np.exp(2j * np.pi * f[None, :] * tau_hat[:, None]) * np.exp(1j * psi_hat)[:, None]

def phase_eq21(hb, f):
    """Baseline (21): LS fit of unwrapped phase.  Model: angle(hb) = angle(b) - 2 pi f tau - psi."""
    P, K = hb.shape
    ph = np.unwrap(np.angle(hb), axis=1)
    A = np.stack([2 * np.pi * f, np.ones(K)], 1)
    sol = np.linalg.lstsq(A, -ph.T, rcond=None)[0]             # [tau; psi] per frame
    return sol[0], sol[1]

def phase_eq22(hb):
    """Baseline (22a),(22b) from IEEE 802.11az."""
    K = hb.shape[1]
    tau = TS / (2 * np.pi) * np.angle(np.sum(hb[:, :-1] * np.conj(hb[:, 1:]), axis=1))
    f = np.arange(K) / TS
    psi = -np.angle(np.sum(hb * np.exp(2j * np.pi * f[None, :] * tau[:, None]), axis=1))
    return tau, psi

def robust_unwrap(ang_k, om_k):
    """Eq.(26): robust unwrap using local sums over k and its 3 preceding / 3 succeeding indices."""
    n = len(ang_k); loc = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - 3), min(n, i + 4)
        loc[i] = np.angle(np.sum(om_k[lo:hi]))
    s = np.unwrap(loc)
    return np.mod(ang_k - s + np.pi, 2 * np.pi) - np.pi + s

def phase_algo4_eq25(hb, f, thr=0.1):
    P, K = hb.shape
    tau_c, psi_c = phase_eq22(hb)                                         # coarse
    bbar = np.mean(hb * np.exp(1j * (2 * np.pi * f[None, :] * tau_c[:, None] + psi_c[:, None])), axis=0)  # (23)
    Kbar = np.where(np.abs(bbar) ** 2 > thr)[0]
    tau_hat = np.empty(P); psi_hat = np.empty(P)
    for p in range(P):
        om = np.conj(hb[p, Kbar]) * bbar[Kbar] * np.exp(-2j * np.pi * f[Kbar] * tau_c[p])
        phi = robust_unwrap(np.angle(om), om)
        w = np.abs(om)
        A = np.stack([2 * np.pi * f[Kbar], np.ones(len(Kbar))], 1)       # unknowns: (tau - tau_c), psi
        Aw = A * np.sqrt(w)[:, None]
        sol = np.linalg.lstsq(Aw, phi * np.sqrt(w), rcond=None)[0]
        tau_hat[p] = tau_c[p] + sol[0]; psi_hat[p] = sol[1]
    return tau_hat, psi_hat

# ---------------------------------------------------------------- metric (Eq. 29)
def chi_metric(hh, sim):
    b, d, f, K = sim["b"], sim["d"], sim["f"], hh.shape[1]
    P = hh.shape[0]; gamma = sim["gamma"]
    bh = hh.mean(0)
    taus = np.linspace(-4e-7, 4e-7, 1601)
    c = np.abs(np.exp(2j * np.pi * f[None, :] * taus[:, None]) @ (b * np.conj(bh)))
    tau = taus[np.argmax(c)]
    num = np.abs(np.sum(np.conj(hh - bh[None, :]) * d * np.exp(2j * np.pi * f[None, :] * tau))) ** 2
    den = (1 - gamma) * K * P * np.sum(np.abs(hh - bh[None, :]) ** 2)
    return float(num / den)

def snr_from_chi(chi):
    return chi / (1 - chi + 1e-12)

# ---------------------------------------------------------------- experiments
def run_phase_trial(sim):
    """ideal gain correction, then compare phase methods"""
    hb = sim["ht"] / sim["g"][:, None]
    f = sim["f"]; out = {}
    t, s = phase_eq21(hb, f);        out["Eq21 (baseline)"] = chi_metric(apply_phase(hb, f, t, s), sim)
    t, s = phase_eq22(hb);           out["Eq22 (802.11az baseline)"] = chi_metric(apply_phase(hb, f, t, s), sim)
    t, s = phase_algo4_eq25(hb, f);  out["Algo4+Eq25 (proposed)"] = chi_metric(apply_phase(hb, f, t, s), sim)
    return out

def run_gain_trial(sim):
    """ideal phase correction, then compare gain methods"""
    f = sim["f"]
    hp = sim["ht"] * np.exp(2j * np.pi * f[None, :] * sim["tau"][:, None]) * np.exp(1j * sim["psi"])[:, None]
    out = {}
    out["Eq6 (baseline)"] = chi_metric(hp / gain_eq6(sim["ht"])[:, None], sim)
    out["Algo1"] = chi_metric(hp / gain_algo1(sim["ht"])[:, None], sim)
    out["Algo3 (proposed)"] = chi_metric(hp / gain_algo3(sim["ht"])[:, None], sim)
    out["No gain correction"] = chi_metric(hp, sim)
    return out

def median_over_trials(fn, n_trials, seed, **kw):
    rng = np.random.default_rng(seed); res = {}
    for _ in range(n_trials):
        sim = simulate(rng=rng, **kw)
        for k, v in fn(sim).items():
            res.setdefault(k, []).append(snr_from_chi(v))
    return {k: (np.median(v), np.percentile(v, 25), np.percentile(v, 75)) for k, v in res.items()}
