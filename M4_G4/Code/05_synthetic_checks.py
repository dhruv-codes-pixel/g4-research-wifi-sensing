"""SYNTHETIC ONLY - known-answer check of csi_methods on the REAL tone map (242 tones, DC gap, Ts = 12.8 us).
Not real-data evidence. Model = paper Eq.(1)-(2): h = b + d (b static with zero-slope phase and mean|b|^2 = gamma; d i.i.d. CN with variance 1-gamma),
observed = g * h * exp(-j2 pi f tau) * exp(-j psi); tau ~ U(-0.55,+0.55) us (covers the 0.2-0.5 us excursions seen in the real data), psi ~ U(-pi,pi),
gain g on two levels 6 dB apart. Known floor of V after perfect correction = (1-gamma)(1-1/P) -> -10.0 dB (gamma 0.9), -20 dB (0.99), -30 dB (0.999).
Usage: python 05_synthetic_checks.py --out ../Results/synthetic"""
import argparse, os, numpy as np, pandas as pd
import csi_methods as M
K, TS, P = M.HE20_KIDX, M.TS_HE, 300
rng = np.random.default_rng(20261114)
def static_channel(kind):
    f = K / TS
    if kind == "flat_LoS":      # strong flat line-of-sight + weak multipath (the paper's Method-1 assumption)
        taps = np.array([1.0, 0.25, 0.15]); dly = np.array([0, 30e-9, 60e-9])
    else:                       # 6-tap exponential multipath, no dominant flat path
        dly = np.arange(6) * 30e-9; taps = np.exp(-dly / 60e-9) ** 0.5
    c = taps * (rng.standard_normal(len(taps)) + 1j * rng.standard_normal(len(taps))) / np.sqrt(2)
    if kind == "flat_LoS": c[0] = taps[0]
    b = (c[None, :] * np.exp(-2j * np.pi * f[:, None] * dly[None, :])).sum(1)
    j = np.where(np.diff(K) == 1)[0]; ts_ = TS / (2 * np.pi) * np.angle(np.sum(b[j] * np.conj(b[j + 1])))
    b = b * np.exp(2j * np.pi * f * ts_)                    # zero-slope condition of the paper
    return b / np.sqrt(np.mean(np.abs(b) ** 2))
rows = []
for kind in ["flat_LoS", "multipath6"]:
    for gamma in [0.9, 0.99, 0.999]:
        res = {m: [] for m in ["A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"]}; terr = {m: [] for m in res}; verr = {m: [] for m in res}
        for trial in range(60):
            b = np.sqrt(gamma) * static_channel(kind)
            d = np.sqrt((1 - gamma) / 2) * (rng.standard_normal((P, len(K))) + 1j * rng.standard_normal((P, len(K))))
            H = b[None, :] + d; tau = rng.uniform(-0.55e-6, 0.55e-6, P); psi = rng.uniform(-np.pi, np.pi, P); g = 10 ** (rng.choice([0, 6.0206], P) / 20)
            X = g[:, None] * H * np.exp(-2j * np.pi * (K / TS)[None, :] * tau[:, None]) * np.exp(-1j * psi)[:, None]
            for m in res:
                Y = M.correct_m4(m, X, K, TS); res[m].append(M.dynamic_power_fraction_db(Y))
        for m in res: rows.append(dict(SYNTHETIC=True, channel=kind, gamma=gamma, method=m, V_dB_mean=np.mean(res[m]), V_dB_min=np.min(res[m]), V_dB_max=np.max(res[m]), known_floor_dB=10 * np.log10((1 - gamma) * (1 - 1 / P)), n_trials=60))
        print(kind, gamma, {m: round(np.mean(v), 2) for m, v in res.items()}, "floor", round(10 * np.log10((1 - gamma) * (1 - 1 / P)), 2), flush=True)
os.makedirs(a_out := (os.sys.argv[os.sys.argv.index("--out") + 1] if "--out" in os.sys.argv else "../Results/synthetic"), exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(a_out, "SYNTHETIC_known_answer_V.csv"), index=False)
