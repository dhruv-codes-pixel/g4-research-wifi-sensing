"""Step 6: figures. REAL_* = real data (Brescia 20 MHz), SYNTHETIC_* = simulation only. Usage: python 06_figures.py --root /path/to/20MHz"""
import argparse, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import csi_methods as M, load_brescia as L
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--run", default="../Results/real_data/run1"); ap.add_argument("--fig", default="../Figures"); a = ap.parse_args()
os.makedirs(a.fig, exist_ok=True); K, TS = M.HE20_KIDX, M.TS_HE
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": .25})
W = pd.read_csv(f"{a.run}/windows.csv"); inv = pd.read_csv("../Data/dataset_inventory.csv"); prof = pd.read_csv("../Data/tone_power_profile.csv")
full = set(inv.loc[inv.full_length == 1, ["cls", "date", "capture"]].apply(tuple, axis=1)); W["full"] = [(c, d, k) in full for c, d, k in zip(W.cls, W.date, W.capture)]
CL = ["E_day", "E_night", "S1", "M4"]; LAB = {"E_day": "E_day\n(0 people, day)", "E_night": "E_night\n(0 people, night)", "S1": "S1\n(1 person sitting)", "M4": "M4\n(4 people, meeting)"}
MS = ["A0", "A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"]; COL = dict(zip(MS, ["#888888", "#d95f02", "#e6ab02", "#1b9e77", "#7570b3"])); NAME = {"A0": "A0 gain only", "A1_eq21": "A1 Eq.(21) as provided", "A1w_eq21_wide": "A1w Eq.(21) wide n", "A2_eq22": "A2 Eq.(22)", "A3_algo4": "A3 Algorithm 4"}

# ---- Fig 1: dataset overview ----
fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.3))
for key, c in [("E_day_capture1_2024-07-12", "#1b9e77"), ("M4_capture0_2024-07-24", "#7570b3")]:
    ax[0].plot(prof.row, prof[key], lw=.9, label=key.split("_capture")[0], color=c)
ax[0].axvspan(-.5, 5.5, color="k", alpha=.12); ax[0].axvspan(126.5, 129.5, color="k", alpha=.12); ax[0].axvspan(250.5, 255.5, color="k", alpha=.12)
ax[0].set_ylim(40, 65); ax[0].set_xlabel("FFT row (shaded = null guard/DC rows, no energy)"); ax[0].set_ylabel("mean power (dB)"); ax[0].set_title("(a) per-tone mean power, 256 rows, 242 usable"); ax[0].legend()
for (cls, p, c) in [("E_day", "0people/2024-07-12/capture1", "#1b9e77"), ("E_night", "0people/2024-07-24/capture1", "#d95f02"), ("M4", "4people/2024-07-24/capture0", "#7570b3")]:
    X = L.load_active(f"{a.root}/{p}.mat")[200:200 + 1500]; fp = 10 * np.log10(np.mean(np.abs(X) ** 2, axis=1)); ax[1].plot(np.arange(len(fp)) * 0.0318, fp - np.median(fp), lw=.6, label=cls, color=c)
ax[1].set_xlabel("time (s)"); ax[1].set_ylabel("frame power - median (dB)"); ax[1].set_title("(b) frame power: discrete 6 dB steps (E_night) before Eq.(6)"); ax[1].legend(loc="upper right")
ii = inv.copy(); ii["lab"] = ii.cls + "\n" + ii.capture.str.replace("capture", "c"); ax[2].bar(range(len(ii)), ii.duration_s / 60, color=["#1b9e77" if f else "#bbbbbb" for f in ii.full_length])
ax[2].set_xticks(range(len(ii))); ax[2].set_xticklabels(ii.lab, rotation=90, fontsize=6); ax[2].set_ylabel("minutes"); ax[2].set_title("(c) 15 captures (grey = partial)")
plt.tight_layout(); plt.savefig(f"{a.fig}/REAL_fig1_dataset_overview.png", dpi=170); plt.close()

# ---- Fig 2: KPI1 per capture ----
cap = W[W.full & W.method.isin(MS)].groupby(["cls", "capture", "method"]).V_dB.mean().reset_index()
fig, ax = plt.subplots(1, 4, figsize=(12.5, 3.4), sharey=True)
for a_, c in zip(ax, CL):
    for j, m in enumerate(MS):
        v = cap[(cap.cls == c) & (cap.method == m)].V_dB.values; a_.scatter(np.full(len(v), j) + np.linspace(-.12, .12, len(v)), v, s=22, color=COL[m], alpha=.85); a_.hlines(v.mean(), j - .3, j + .3, color=COL[m], lw=2)
    a_.set_xticks(range(len(MS))); a_.set_xticklabels(["A0", "A1", "A1w", "A2", "A3"]); a_.set_title(LAB[c].replace("\n", " "), fontsize=8)
ax[0].set_ylabel("KPI1  V (dB), lower = better\n(dot = capture, bar = mean)"); plt.tight_layout(); plt.savefig(f"{a.fig}/REAL_fig2_KPI1_V_by_class.png", dpi=170); plt.close()

# ---- Fig 3: phase of one subcarrier vs packet ----
def pick(path, w):
    X = L.load_active(f"{a.root}/{path}.mat")[200:][w * 300:(w + 1) * 300]; return X
fig, axs = plt.subplots(2, 2, figsize=(12.5, 5.6))
for r, (title, path, w) in enumerate([("E_day capture1, window 6 (contains Eq.21 clipped-offset frames)", "0people/2024-07-12/capture1", 6), ("M4 capture0, window 10 (4 people)", "4people/2024-07-24/capture0", 10)]):
    Xw = pick(path, w); kk = 60
    for m in ["A0", "A1_eq21", "A2_eq22", "A3_algo4"]:
        Y = M.correct_m4(m, Xw, K, TS); mu = Y[:, kk].mean(); ph = np.angle(Y[:, kk] * np.conj(mu))
        axs[r, 0 if m == "A0" else 1].plot(ph, lw=.8, color=COL[m], label=NAME[m]) if m != "A0" else axs[r, 0].plot(ph, lw=.5, color=COL[m], label=NAME[m])
    axs[r, 0].set_title(f"{title}\nA0 gain only: phase of tone k={K[kk]} minus its mean", fontsize=8); axs[r, 1].set_title("after phase correction (same tone)", fontsize=8); axs[r, 1].legend(fontsize=7, ncol=3, loc="upper right")
    for c_ in (0, 1): axs[r, c_].set_xlabel("packet index in window"); axs[r, c_].set_ylabel("phase (rad)")
    axs[r, 1].set_ylim(-0.6, 0.6)
plt.tight_layout(); plt.savefig(f"{a.fig}/REAL_fig3_phase_vs_packet.png", dpi=170); plt.close()

# ---- Fig 4: Eq.(21) clipped offset search ----
X = L.load_active(f"{a.root}/0people/2024-07-12/capture1.mat")[200:]; xs, y1, y2 = [], [], []
for w in range(0, 40):
    G = M.gain_eq6(X[w * 300:(w + 1) * 300]); t1, _ = M.eq21(G, K, TS); t1w, _ = M.eq21_wide(G, K, TS); t3, _ = M.algo4_eq25(G, K, TS)
    xs += list((t3 - np.median(t3)) * 1e9); y1 += list(((t1 - np.median(t1)) - (t3 - np.median(t3))) * 1e9); y2 += list(((t1w - np.median(t1w)) - (t3 - np.median(t3))) * 1e9)
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.3), sharey=True)
ax[0].scatter(xs, y1, s=3, color=COL["A1_eq21"]); ax[0].set_title("A1: Eq.(21) as provided (offset n searched in -3..3)", fontsize=8.5); ax[1].scatter(xs, y2, s=3, color=COL["A1w_eq21_wide"]); ax[1].set_title("A1w: same fit, n searched in -8..8", fontsize=8.5)
for a_ in ax: a_.set_xlabel("timing excursion of the frame, tau(A3) - window median (ns)"); a_.axhline(0, color="k", lw=.5)
ax[0].set_ylabel("tau(Eq.21) - tau(Algorithm 4) (ns)"); plt.suptitle("E_day capture1, 40 windows: errors appear only on frames with large timing excursions and come in steps of about 78 ns", fontsize=9, y=1.02)
plt.tight_layout(); plt.savefig(f"{a.fig}/REAL_fig4_eq21_offset_search.png", dpi=170, bbox_inches="tight"); plt.close()

# ---- Fig 5: window-level V of A3 per capture ----
fig, ax = plt.subplots(1, 2, figsize=(11.5, 3.4)); wv = W[W.full & (W.method == "A3_algo4")]
pos = 0; ticks = []; labs = []
for c in CL:
    for cp in sorted(wv[wv.cls == c].capture.unique()):
        v = wv[(wv.cls == c) & (wv.capture == cp)].V_dB.values; ax[0].boxplot(v, positions=[pos], widths=.6, showfliers=False, patch_artist=True, boxprops=dict(facecolor={"E_day": "#1b9e77", "E_night": "#d95f02", "S1": "#e6ab02", "M4": "#7570b3"}[c], alpha=.6)); ticks.append(pos); labs.append(cp[-1]); pos += 1
    pos += .7
ax[0].set_xticks(ticks); ax[0].set_xticklabels(labs); ax[0].set_ylabel("A3 window V (dB)"); ax[0].set_title("A3: window-level V per capture (colours = class: E_day, E_night, S1, M4)", fontsize=8.5); ax[0].set_xlabel("capture index")
tm = W[W.method.isin(MS)].groupby("method").time_ms.mean().reindex(MS); ax[1].bar(range(len(MS)), tm.values, color=[COL[m] for m in MS]); ax[1].set_xticks(range(len(MS))); ax[1].set_xticklabels(["A0", "A1", "A1w", "A2", "A3"]); ax[1].set_ylabel("ms per 300-frame window (1 CPU, incl. Eq.6)"); ax[1].set_title("KPI3: run time", fontsize=8.5)
plt.tight_layout(); plt.savefig(f"{a.fig}/REAL_fig5_window_V_and_time.png", dpi=170); plt.close()

# ---- Fig 6: synthetic ----
S = pd.read_csv("../Results/synthetic/SYNTHETIC_known_answer_V.csv"); fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.3), sharey=True)
for a_, kind in zip(ax, ["flat_LoS", "multipath6"]):
    x = np.arange(3)
    for j, m in enumerate(["A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"]):
        v = S[(S.channel == kind) & (S.method == m)].sort_values("gamma"); a_.plot(x, v.V_dB_mean, "o-", color=COL[m], label=NAME[m])
    fl = S[(S.channel == kind) & (S.method == "A3_algo4")].sort_values("gamma").known_floor_dB; a_.plot(x, fl, "k--", label="known floor 10log10(1-gamma)"); a_.set_xticks(x); a_.set_xticklabels(["0.9", "0.99", "0.999"]); a_.set_xlabel("static fraction gamma"); a_.set_title(f"SYNTHETIC ONLY: {kind}", fontsize=8.5)
ax[0].set_ylabel("V (dB)"); ax[0].legend(fontsize=7); plt.tight_layout(); plt.savefig(f"{a.fig}/SYNTHETIC_fig6_known_answer_V.png", dpi=170); plt.close()
print(sorted(os.listdir(a.fig)))
