"""
M3 Member E: controlled single-reflector Doppler experiment.

Keep this file beside csi_sim.py.

Quick smoke test:
    python doppler_exp.py 2

Full assigned experiment:
    python doppler_exp.py 200

The simulation uses a direct 20 Hz complex sinusoid as a controlled
single-reflector model. It is not real hardware CSI.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import csi_sim as c

DOPPLER_HZ = 20.0
FRAME_INTERVAL = 0.005
N_FRAMES = 512
GAMMA = 0.9
K_VALUES = [64, 256, 512]
TIMING_MAX_VALUES = [100e-9, 800e-9]
DEFAULT_TRIALS = 200
BASE_SEED = 2026
NFFT = 4096

METHODS = [
    "No correction",
    "Eq6 + Eq22",
    "Eq6 + Eq21",
    "Eq6 + Algorithm4",
    "Ideal correction",
]

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def make_doppler_trial(K, tau_max, rng):
    f = np.arange(K) / c.TS
    t = np.arange(N_FRAMES) * FRAME_INTERVAL

    b = c.make_static(K, GAMMA, rng)

    # One reflector: one temporal complex tone at 20 Hz.
    reflector_delay = 80e-9
    reflector_response = np.exp(-2j * np.pi * f * reflector_delay)
    doppler_wave = np.exp(2j * np.pi * DOPPLER_HZ * t)
    d = (
        np.sqrt(1.0 - GAMMA)
        * doppler_wave[:, None]
        * reflector_response[None, :]
    )
    h = b[None, :] + d

    # Receiver gain: slow fluctuation plus discrete AGC.
    w = rng.standard_normal(N_FRAMES + 120)
    kernel = np.ones(21) / 21
    g1 = np.convolve(w, kernel, mode="same")[60:60 + N_FRAMES]
    if g1.std() > 0:
        g1 = 0.2 * g1 / g1.std()
    g2 = rng.choice([-0.5, 0.0, 0.5], size=N_FRAMES, p=[0.2, 0.6, 0.2])
    gain = 10 ** ((g1 + g2) / 20.0)

    tau = rng.uniform(0.0, tau_max, N_FRAMES)
    psi = rng.uniform(-np.pi, np.pi, N_FRAMES)

    ht = (
        gain[:, None] * h
        * np.exp(-2j * np.pi * f[None, :] * tau[:, None])
        * np.exp(-1j * psi)[:, None]
    )

    return {
        "ht": ht, "h": h, "b": b, "d": d,
        "g": gain, "tau": tau, "psi": psi,
        "f": f, "t": t, "reflector_response": reflector_response,
        "gamma": GAMMA, "K": K,
    }


def apply_methods(sim):
    ht, f = sim["ht"], sim["f"]
    corrected = {"No correction": ht}

    # Eq. 6 + Eq. 22
    ghat = c.gain_eq6(ht)
    hb = ht / ghat[:, None]
    tau_hat, psi_hat = c.phase_eq22(hb)
    corrected["Eq6 + Eq22"] = c.apply_phase(hb, f, tau_hat, psi_hat)

    # Eq. 6 + Eq. 21
    tau_hat, psi_hat = c.phase_eq21(hb, f)
    corrected["Eq6 + Eq21"] = c.apply_phase(hb, f, tau_hat, psi_hat)

    # Eq. 6 + Algorithm 4 + Eq. 25
    tau_hat, psi_hat = c.phase_algo4_eq25(hb, f)
    corrected["Eq6 + Algorithm4"] = c.apply_phase(hb, f, tau_hat, psi_hat)

    # Ideal receiver-error correction
    corrected["Ideal correction"] = (
        ht / sim["g"][:, None]
        * np.exp(2j * np.pi * f[None, :] * sim["tau"][:, None])
        * np.exp(1j * sim["psi"])[:, None]
    )
    return corrected


def _dominant_temporal_component(dynamic):
    """
    Estimate the dominant temporal component from CSI alone.

    Uses a randomized rank-1 SVD approximation. It does not use the
    simulated reflector frequency response or any ground-truth channel
    vector. The singular-vector phase ambiguity does not affect the
    magnitude spectrum.
    """
    x = np.asarray(dynamic, dtype=np.complex128)
    n_frames, n_subcarriers = x.shape

    # Randomized range finder with a small oversampling dimension.
    rank = min(8, n_frames, n_subcarriers)
    rng = np.random.default_rng(731)
    omega = (
        rng.standard_normal((n_subcarriers, rank))
        + 1j * rng.standard_normal((n_subcarriers, rank))
    ) / np.sqrt(2.0)

    y = x @ omega
    q, _ = np.linalg.qr(y, mode="reduced")

    # Two power iterations improve approximation of the leading component.
    for _ in range(2):
        y = x @ (x.conj().T @ q)
        q, _ = np.linalg.qr(y, mode="reduced")

    small = q.conj().T @ x
    u_small, singular_values, _ = np.linalg.svd(
        small, full_matrices=False
    )
    temporal = q @ u_small[:, 0]
    return temporal * singular_values[0]


def doppler_metrics(hh, sim=None):
    """
    Estimate Doppler metrics from corrected CSI without oracle combining.

    The static component is removed by temporal centering. A rank-1
    temporal component is then estimated from the CSI matrix itself using
    randomized SVD. The known reflector_response in `sim` is deliberately
    not used.
    """
    dynamic = hh - np.mean(hh, axis=0, keepdims=True)
    temporal = _dominant_temporal_component(dynamic)
    temporal = temporal - np.mean(temporal)

    window = np.hanning(len(temporal))
    spec = np.fft.fft(temporal * window, n=NFFT)
    freqs = np.fft.fftfreq(NFFT, d=FRAME_INTERVAL)
    power = np.abs(spec) ** 2

    positive = freqs >= 0
    fpos, ppos = freqs[positive], power[positive]
    peak_i = int(np.argmax(ppos))
    peak_freq = float(fpos[peak_i])
    peak_error = abs(peak_freq - DOPPLER_HZ)

    # Median floor excluding bins within +/- 3 bins of the detected peak.
    mask = np.ones(len(ppos), dtype=bool)
    mask[max(0, peak_i - 3):min(len(ppos), peak_i + 4)] = False
    floor = float(np.median(ppos[mask])) if np.any(mask) else 1e-15
    peak_power = float(ppos[peak_i])
    lfr_db = 10 * np.log10((peak_power + 1e-15) / (floor + 1e-15))

    return peak_error, lfr_db, fpos, ppos


def run_experiment(n_trials=DEFAULT_TRIALS):
    if n_trials < 2:
        raise ValueError("Use at least 2 trials.")

    configurations = [
        (K, tau_max)
        for K in K_VALUES
        for tau_max in TIMING_MAX_VALUES
    ]
    peak_errors = np.zeros((len(configurations), n_trials, len(METHODS)))
    lfr_values = np.zeros_like(peak_errors)
    spectra = {}

    for ci, (K, tau_max) in enumerate(configurations):
        print(
            f"\nConfiguration {ci + 1}/{len(configurations)}: "
            f"K={K}, max timing error={tau_max * 1e9:.0f} ns",
            flush=True,
        )

        for trial in range(n_trials):
            seed = BASE_SEED + ci * 100000 + trial
            rng = np.random.default_rng(seed)
            sim = make_doppler_trial(K, tau_max, rng)
            corrected = apply_methods(sim)

            for mi, method in enumerate(METHODS):
                err, lfr, fpos, ppos = doppler_metrics(corrected[method], sim)
                peak_errors[ci, trial, mi] = err
                lfr_values[ci, trial, mi] = lfr
                if trial == 0:
                    spectra[(ci, method)] = (fpos.copy(), ppos.copy())

            if (trial + 1) % 10 == 0 or trial == 0:
                print(f"  completed {trial + 1}/{n_trials} trials", flush=True)

    np.savez_compressed(
        OUT_DIR / "doppler_results.npz",
        configurations=np.asarray(configurations, dtype=float),
        methods=np.asarray(METHODS),
        peak_errors_hz=peak_errors,
        line_to_floor_db=lfr_values,
        n_trials=n_trials,
        doppler_hz=DOPPLER_HZ,
        frame_interval_s=FRAME_INTERVAL,
        n_frames=N_FRAMES,
        gamma=GAMMA,
        seed=BASE_SEED,
    )

    write_summary(configurations, peak_errors, lfr_values)
    plot_summary(configurations, peak_errors, lfr_values)
    plot_spectra(configurations, spectra)

    print("\nExperiment complete.", flush=True)
    print(f"Results saved in: {OUT_DIR}", flush=True)


def write_summary(configurations, peak_errors, lfr_values):
    path = OUT_DIR / "doppler_summary.csv"
    with path.open("w", encoding="utf-8") as file:
        file.write(
            "K,timing_max_ns,method,peak_error_median_hz,"
            "peak_error_q2.5_hz,peak_error_q97.5_hz,"
            "lfr_median_db,lfr_q2.5_db,lfr_q97.5_db\n"
        )
        for ci, (K, tau_max) in enumerate(configurations):
            for mi, method in enumerate(METHODS):
                e = peak_errors[ci, :, mi]
                l = lfr_values[ci, :, mi]
                file.write(
                    f"{K},{tau_max * 1e9:.0f},{method},"
                    f"{np.median(e):.6f},{np.percentile(e, 2.5):.6f},"
                    f"{np.percentile(e, 97.5):.6f},"
                    f"{np.median(l):.6f},{np.percentile(l, 2.5):.6f},"
                    f"{np.percentile(l, 97.5):.6f}\n"
                )
    print(f"Saved summary: {path}", flush=True)


def plot_summary(configurations, peak_errors, lfr_values):
    labels = [f"K={K}\n{tau * 1e9:.0f} ns" for K, tau in configurations]
    x = np.arange(len(configurations))
    width = 0.14

    plt.figure(figsize=(12, 5))
    for mi, method in enumerate(METHODS):
        med = np.median(peak_errors[:, :, mi], axis=1)
        plt.bar(x + (mi - 2) * width, med, width, label=method)
    plt.xticks(x, labels)
    plt.ylabel("Median peak-frequency error (Hz)")
    plt.xlabel("Configuration")
    plt.title("Doppler peak-frequency error comparison")
    plt.legend(fontsize=8)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "doppler_peak_error.png", dpi=200)
    plt.close()

    plt.figure(figsize=(12, 5))
    for mi, method in enumerate(METHODS):
        med = np.median(lfr_values[:, :, mi], axis=1)
        plt.bar(x + (mi - 2) * width, med, width, label=method)
    plt.xticks(x, labels)
    plt.ylabel("Median line-to-floor ratio (dB)")
    plt.xlabel("Configuration")
    plt.title("Doppler line-to-floor ratio comparison")
    plt.legend(fontsize=8)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "doppler_line_to_floor.png", dpi=200)
    plt.close()


def plot_spectra(configurations, spectra):
    for ci, (K, tau_max) in enumerate(configurations):
        plt.figure(figsize=(9, 5))
        for method in METHODS:
            freqs, power = spectra[(ci, method)]
            power_db = 10 * np.log10(
                power / (np.max(power) + 1e-15) + 1e-15
            )
            plt.plot(freqs, power_db, label=method)
        plt.axvline(20.0, linestyle="--", label="True Doppler = 20 Hz")
        plt.xlim(0, 50)
        plt.ylim(-80, 5)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Normalised spectrum (dB)")
        plt.title(
            f"Doppler spectrum: K={K}, timing max={tau_max * 1e9:.0f} ns"
        )
        plt.grid(alpha=0.3)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(
            OUT_DIR / f"doppler_spectrum_K{K}_{int(tau_max * 1e9)}ns.png",
            dpi=200,
        )
        plt.close()


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_TRIALS
    run_experiment(n)
