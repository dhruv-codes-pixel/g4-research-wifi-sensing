"""
M2 CSI simulator and correction algorithms for the M3 Doppler experiment.

Keep this file in the same folder as doppler_exp.py.
This is the complete simulation/correction module used by that experiment.
"""

import numpy as np
from sklearn.cluster import DBSCAN
from scipy.stats import norm

# OFDM symbol duration and frame repetition interval
TS = 3.2e-6
TREP = 0.1


# ---------------------------------------------------------------------
# Static channel and CSI simulator
# ---------------------------------------------------------------------
def make_static(K, gamma, rng, max_delay=200e-9, n_taps=6):
    """Generate a static multipath channel with normalized power gamma."""
    f = np.arange(K) / TS
    delays = np.sort(rng.uniform(0, max_delay, n_taps))
    delays[0] = 0.0

    amps = (
        rng.standard_normal(n_taps)
        + 1j * rng.standard_normal(n_taps)
    )
    amps *= np.exp(-delays / 80e-9)
    amps[0] *= 3.0

    b = (
        amps[None, :]
        * np.exp(-2j * np.pi * f[:, None] * delays[None, :])
    ).sum(axis=1)

    # Enforce the zero-slope phase constraint.
    theta = np.angle(np.sum(b[:-1] * np.conj(b[1:])))
    b = b * np.exp(1j * theta * np.arange(K))

    # Normalize average static-channel power to gamma.
    b *= np.sqrt(gamma * K / np.sum(np.abs(b) ** 2))
    return b


def simulate(
    K=64,
    P=300,
    gamma=0.9,
    tau_max=1e-7,
    rng=None,
    agc=(-0.5, 0.0, 0.5),
    agc_p=(0.2, 0.6, 0.2),
):
    """Simulate CSI with dynamic component, gain, timing and phase errors."""
    if rng is None:
        rng = np.random.default_rng()

    f = np.arange(K) / TS
    b = make_static(K, gamma, rng)

    d = np.sqrt((1 - gamma) / 2) * (
        rng.standard_normal((P, K))
        + 1j * rng.standard_normal((P, K))
    )
    h = b[None, :] + d

    # Slow gain variation and discrete AGC steps.
    w = rng.standard_normal(P + 120)
    kernel = np.ones(21) / 21
    g1 = np.convolve(w, kernel, mode="same")[60:60 + P]
    if g1.std() > 0:
        g1 = 0.2 * g1 / g1.std()

    g2 = rng.choice(agc, size=P, p=agc_p)
    g = 10 ** ((g1 + g2) / 20)

    tau = rng.uniform(0, tau_max, P)
    psi = rng.uniform(-np.pi, np.pi, P)

    ht = (
        g[:, None]
        * h
        * np.exp(-2j * np.pi * f[None, :] * tau[:, None])
        * np.exp(-1j * psi)[:, None]
    )

    return {
        "ht": ht,
        "h": h,
        "b": b,
        "d": d,
        "g": g,
        "tau": tau,
        "psi": psi,
        "f": f,
        "gamma": gamma,
    }


# ---------------------------------------------------------------------
# Gain correction
# ---------------------------------------------------------------------
def gain_eq6(ht):
    """Baseline gain estimate (Equation 6)."""
    return np.sqrt(np.mean(np.abs(ht) ** 2, axis=1))


def lpf_ma(x, half_width):
    """Moving-average low-pass filter."""
    n = len(x)
    out = np.empty(n)

    for i in range(n):
        lo = max(0, i - half_width)
        hi = min(n, i + half_width + 1)
        out[i] = x[lo:hi].mean()

    return out


def gain_algo1(ht, eps=0.2):
    """Gain correction Algorithm 1: DBSCAN plus low-pass filtering."""
    P = ht.shape[0]
    G = 10 * np.log10(np.mean(np.abs(ht) ** 2, axis=1))

    dG = np.zeros(P)
    dG[1:] = G[1:] - G[:-1]

    labels = DBSCAN(eps=eps, min_samples=1).fit_predict(
        dG[1:].reshape(-1, 1)
    )
    cluster_means = {
        label: dG[1:][labels == label].mean()
        for label in np.unique(labels)
    }

    g2 = np.zeros(P)
    for p in range(1, P):
        g2[p] = g2[p - 1] + cluster_means[labels[p - 1]]

    g1 = lpf_ma(G - g2, half_width=int(round(6 / TREP)))
    return 10 ** ((g1 + g2) / 20)


def _D(x, zmax=30):
    """Auxiliary objective function used by Algorithm 3."""
    z = np.arange(-zmax, zmax + 1)
    return np.sum(
        (
            norm.sf((z - 0.5) * x)
            - norm.sf((z + 0.5) * x)
        ) * z ** 2
    )


def gain_algo3(ht, hw=None):
    """Algorithm 3: ML gain estimation for a uniformly spaced AGC grid."""
    if hw is None:
        hw = int(round(6 / TREP))

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
        sig2 = -(lam ** 2) / (2 * np.pi ** 2) * np.log(
            max(m, 1e-12)
        )
        objective = (
            sig2
            + lam ** 2 * _D(
                lam / np.sqrt(max(sig2, 1e-12))
            )
        )

        if objective < best[0]:
            best = (objective, 10 ** ((g1 + g2) / 20))

    return best[1] if best[1] is not None else gain_eq6(ht)


# ---------------------------------------------------------------------
# Phase correction
# ---------------------------------------------------------------------
def apply_phase(hb, f, tau_hat, psi_hat):
    """Apply estimated timing and common phase correction."""
    return (
        hb
        * np.exp(2j * np.pi * f[None, :] * tau_hat[:, None])
        * np.exp(1j * psi_hat)[:, None]
    )


def phase_eq21(hb, f):
    """Baseline Equation 21: least-squares phase-slope estimation."""
    P, K = hb.shape
    ph = np.unwrap(np.angle(hb), axis=1)

    A = np.stack([2 * np.pi * f, np.ones(K)], axis=1)
    solution = np.linalg.lstsq(A, -ph.T, rcond=None)[0]

    return solution[0], solution[1]


def phase_eq22(hb):
    """Baseline Equations 22a and 22b."""
    K = hb.shape[1]

    tau = (
        TS / (2 * np.pi)
        * np.angle(
            np.sum(
                hb[:, :-1] * np.conj(hb[:, 1:]),
                axis=1,
            )
        )
    )

    f = np.arange(K) / TS

    psi = -np.angle(
        np.sum(
            hb * np.exp(
                2j * np.pi * f[None, :] * tau[:, None]
            ),
            axis=1,
        )
    )

    return tau, psi


def robust_unwrap(ang_k, om_k):
    """Robust phase unwrapping using local sums (Equation 26)."""
    n = len(ang_k)
    local_phase = np.empty(n)

    for i in range(n):
        lo = max(0, i - 3)
        hi = min(n, i + 4)
        local_phase[i] = np.angle(np.sum(om_k[lo:hi]))

    smooth_phase = np.unwrap(local_phase)

    return (
        np.mod(
            ang_k - smooth_phase + np.pi,
            2 * np.pi,
        )
        - np.pi
        + smooth_phase
    )


def phase_algo4_eq25(hb, f, thr=0.1):
    """Proposed Algorithm 4 with Equation 25 refinement."""
    P, K = hb.shape

    tau_coarse, psi_coarse = phase_eq22(hb)

    bbar = np.mean(
        hb
        * np.exp(
            1j * (
                2 * np.pi * f[None, :] * tau_coarse[:, None]
                + psi_coarse[:, None]
            )
        ),
        axis=0,
    )

    selected = np.where(np.abs(bbar) ** 2 > thr)[0]

    # Fallback if too few usable subcarriers remain.
    if len(selected) < 2:
        return tau_coarse, psi_coarse

    tau_hat = np.empty(P)
    psi_hat = np.empty(P)

    for p in range(P):
        om = (
            np.conj(hb[p, selected])
            * bbar[selected]
            * np.exp(-2j * np.pi * f[selected] * tau_coarse[p])
        )

        phi = robust_unwrap(np.angle(om), om)
        weights = np.abs(om)

        A = np.stack(
            [
                2 * np.pi * f[selected],
                np.ones(len(selected)),
            ],
            axis=1,
        )

        Aw = A * np.sqrt(weights)[:, None]
        solution = np.linalg.lstsq(
            Aw,
            phi * np.sqrt(weights),
            rcond=None,
        )[0]

        tau_hat[p] = tau_coarse[p] + solution[0]
        psi_hat[p] = solution[1]

    return tau_hat, psi_hat
