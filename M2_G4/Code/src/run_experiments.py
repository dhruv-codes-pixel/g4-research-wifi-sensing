"""Run M2 experiments. Usage: python run_experiments.py [n_trials]   (default 40; paper used 2000)"""
import sys, json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csi_sim as c

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
OUT = "../../Results/"
res = {}

def curve(fn, key_name, values, **fixed):
    table = {}
    for v in values:
        kw = dict(fixed); kw[key_name] = v
        r = c.median_over_trials(fn, N, seed=123, **kw)
        for m, (med, q1, q3) in r.items():
            table.setdefault(m, []).append((v, med, q1, q3))
    return table

def plot(table, xlabel, title, fname, xscale=None, db=True):
    plt.figure(figsize=(6, 4))
    for m, rows in table.items():
        x = [r[0] for r in rows]; y = [10 * np.log10(r[1]) if db else r[1] for r in rows]
        lo = [10 * np.log10(r[2]) if db else r[2] for r in rows]; hi = [10 * np.log10(r[3]) if db else r[3] for r in rows]
        plt.plot(x, y, marker="o", label=m); plt.fill_between(x, lo, hi, alpha=0.15)
    plt.xlabel(xlabel); plt.ylabel("Median SNR of recovered dynamic CSI (dB)"); plt.title(title)
    if xscale: plt.xscale(xscale)
    plt.grid(alpha=.3); plt.legend(fontsize=7); plt.tight_layout(); plt.savefig(OUT + fname, dpi=150); plt.close()

# E1: reproduction -- phase methods vs gamma (cf. paper Fig. 7)
gam = [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
t = curve(c.run_phase_trial, "gamma", gam, K=64); res["E1_phase_vs_gamma"] = t
plot(t, "gamma (static-channel power fraction)", "E1: phase correction vs gamma (K=64)", "E1_phase_vs_gamma.png")
# E2: reproduction -- gain methods vs gamma (cf. paper Fig. 6)
t = curve(c.run_gain_trial, "gamma", gam, K=64); res["E2_gain_vs_gamma"] = t
plot(t, "gamma", "E2: gain correction vs gamma (K=64)", "E2_gain_vs_gamma.png")
# E3: extension -- timing-error range (stress test beyond paper's 100 ns)
tm = [100e-9, 400e-9, 800e-9, 1.6e-6, 2.4e-6, 3.2e-6]
table = {}
for tmax in tm:
    r = c.median_over_trials(c.run_phase_trial, N, seed=123, K=64, gamma=0.9, tau_max=tmax)
    for m, (med, q1, q3) in r.items(): table.setdefault(m, []).append((tmax * 1e9, med, q1, q3))
res["E3_phase_vs_timing_range_ns"] = table
plot(table, "max timing error tau_max (ns)", "E3: timing-error stress test (K=64, gamma=0.9)", "E3_phase_vs_tau_max.png", xscale="log")
# E4: extension -- bandwidth scaling (number of subcarriers at fixed 312.5 kHz spacing)
t = curve(c.run_phase_trial, "K", [64, 128, 256, 512], gamma=0.9); res["E4_phase_vs_K"] = t
plot(t, "subcarriers K (64=20 MHz, 128=40, 256=80, 512=160 MHz)", "E4: bandwidth scaling, phase correction", "E4_phase_vs_K.png", xscale="log")
t = curve(c.run_gain_trial, "K", [64, 128, 256, 512], gamma=0.9); res["E5_gain_vs_K"] = t
plot(t, "subcarriers K", "E5: bandwidth scaling, gain correction", "E5_gain_vs_K.png", xscale="log")

json.dump(res, open(OUT + "results.json", "w"), indent=1, default=float)
# print summary table
for k, tb in res.items():
    print("\n" + k)
    for m, rows in tb.items():
        print(f"  {m:28s}", " ".join(f"{r[0]:g}:{10*np.log10(r[1]):.1f}dB" for r in rows))
