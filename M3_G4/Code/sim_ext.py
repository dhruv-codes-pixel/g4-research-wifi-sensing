"""
sim_ext.py  --  M3 simulator extensions (csi_sim.py is NOT modified).

simulate(...) reproduces csi_sim.simulate() draw-for-draw when g1_mode="m2" and dyn="i"
(bit-identical output for the same Generator state; checked in tests/test_gain.py), and adds:

  g1_mode="paper" : large-scale gain g1 band-limited to [0, 0.1] Hz as stated in Ratnam et al.
                    Sec. V-B ("Doppler power spectrum with non-zero support on [0, 0.1] Hz").
                    csi_sim.py uses a 21-tap moving average at Trep = 0.1 s, i.e. a 2.1 s window whose
                    first spectral null is at ~0.48 Hz, so its g1 contains power far above 0.1 Hz.
                    Std is normalised to 0.2 dB per realisation exactly as in csi_sim.py, so the ONLY
                    difference from "m2" is the spectral support.
  g1_mode="none"  : g1 = 0 (diagnostic only).
  dyn="ii"        : Doppler-structured dynamic component of Sec. V-B type (ii):
                    d[p,k] = alpha_p * exp(-j 2 pi f_k (tau + tau0)), tau ~ Uni[0, 300 ns),
                    tau0 = LoS delay of b_k, alpha_p complex Gaussian process of variance 1-gamma with
                    Doppler support [0.5, 1] Hz.

Extra random draws for "paper"/"ii" are taken AFTER all csi_sim draws, so b, g2, tau_p, psi_p (and d for
type i) are identical across modes for the same per-trial generator -> paired comparisons.
"""
import numpy as np
import csi_sim as c

TS, TREP = c.TS, c.TREP


def _bandlimited(n, rng, f_lo, f_hi, trep=TREP, n_long=4096, complex_=False):
    """Gaussian process with flat PSD on [f_lo, f_hi] Hz (positive freqs only if complex_)."""
    if complex_:
        w = rng.standard_normal(n_long) + 1j * rng.standard_normal(n_long)
        W = np.fft.fft(w); f = np.fft.fftfreq(n_long, trep)
        W[(f < f_lo) | (f > f_hi)] = 0
        x = np.fft.ifft(W)
    else:
        w = rng.standard_normal(n_long)
        W = np.fft.rfft(w); f = np.fft.rfftfreq(n_long, trep)
        W[(f < f_lo) | (f > f_hi)] = 0
        x = np.fft.irfft(W, n_long)
    s = rng.integers(0, n_long - n)
    return x[s:s + n]


def los_delay(b, f, span=TS / 2, n=4001):
    taus = np.linspace(-span, span, n)
    return taus[np.argmax(np.abs(np.exp(2j * np.pi * f[None, :] * taus[:, None]) @ b))]


def simulate(K=64, P=300, gamma=0.9, tau_max=1e-7, rng=None, agc=(-0.5, 0, 0.5), agc_p=(0.2, 0.6, 0.2),
             g1_mode="m2", dyn="i"):
    rng = rng or np.random.default_rng()
    f = np.arange(K) / TS
    b = c.make_static(K, gamma, rng)
    d = np.sqrt((1 - gamma) / 2) * (rng.standard_normal((P, K)) + 1j * rng.standard_normal((P, K)))
    w = rng.standard_normal(P + 120); ker = np.ones(21) / 21
    g1 = np.convolve(w, ker, mode="same")[60:60 + P]; g1 = 0.2 * g1 / g1.std()
    g2 = rng.choice(agc, size=P, p=agc_p)
    tau = rng.uniform(0, tau_max, P)
    psi = rng.uniform(-np.pi, np.pi, P)
    # ---- extensions (drawn after every csi_sim draw) ----
    if g1_mode == "paper":
        g1 = _bandlimited(P, rng, 0.0, 0.1); g1 = 0.2 * g1 / g1.std()
    elif g1_mode == "none":
        g1 = np.zeros(P)
    elif g1_mode != "m2":
        raise ValueError(g1_mode)
    if dyn == "ii":
        tau_m = rng.uniform(0, 300e-9); tau0 = los_delay(b, f)
        alpha = _bandlimited(P, rng, 0.5, 1.0, complex_=True)
        alpha *= np.sqrt((1 - gamma) / np.mean(np.abs(alpha) ** 2))
        d = alpha[:, None] * np.exp(-2j * np.pi * f * (tau_m + tau0))[None, :]
    elif dyn != "i":
        raise ValueError(dyn)
    h = b[None, :] + d
    g = 10 ** ((g1 + g2) / 20)
    ht = g[:, None] * h * np.exp(-2j * np.pi * f[None, :] * tau[:, None]) * np.exp(-1j * psi)[:, None]
    return dict(ht=ht, h=h, b=b, d=d, g=g, g1=g1, g2=g2, tau=tau, psi=psi, f=f, gamma=gamma)
