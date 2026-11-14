"""Step 1: complete inventory of the Brescia 20 MHz dataset -> Data/dataset_inventory.csv, Data/tone_power_profile.csv.
Usage: python 01_inventory.py --root /path/to/20MHz --out ../Data"""
import argparse, os, numpy as np, pandas as pd
import csi_methods as M, load_brescia as L

ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", default="../Data"); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
rows, prof = [], {}
for c in L.list_captures(a.root):
    Xf = L.load_mat(c["mat"]); ts = L.pcap_timestamps(c["pcap"]); n = Xf.shape[0]
    pw_tone = np.mean(np.abs(Xf) ** 2, axis=0)
    X = Xf[:, M.HE20_ROWS]; G = M.gain_eq6(X)
    dt = np.diff(ts); fp = 10 * np.log10(np.mean(np.abs(X) ** 2, axis=1))
    tau, psi = M.eq22(G, M.HE20_KIDX)
    # raw-vs-processed evidence: sanitised CSI would have psi ~ constant and tau ~ 0; uniform psi => random common phase per frame
    R_psi = abs(np.mean(np.exp(1j * psi)))
    d = np.diff(fp)
    gaps = np.where(dt > 0.5)[0]
    # usable windows: skip first 200 frames, P=300 non-overlapping, no gap > 0.5 s inside
    good = 0
    for w in range((n - 200) // 300):
        s, e = 200 + w * 300, 200 + (w + 1) * 300
        if np.all(dt[s:e - 1] <= 0.5): good += 1
    rows.append(dict(cls=c["cls"], people=c["people"], date=c["date"], capture=c["capture"], code=c["code"], frames_mat=n, packets_pcap=len(ts),
        duration_s=round(ts[-1] - ts[0], 1), dt_median_ms=round(np.median(dt) * 1e3, 2), dt_p1_ms=round(np.percentile(dt, 1) * 1e3, 2), dt_p99_ms=round(np.percentile(dt, 99) * 1e3, 2),
        n_gaps_gt_0p5s=len(gaps), max_gap_s=round(dt.max(), 2), full_length=int(n >= 17000), usable_windows_P300=good,
        frame_power_median_dB=round(np.median(fp), 2), frame_power_std_dB=round(fp.std(), 3), absdiff_p50_dB=round(np.percentile(np.abs(d), 50), 3),
        absdiff_p99_dB=round(np.percentile(np.abs(d), 99), 3), frames_gt3dB_from_median_pct=round(100 * np.mean(np.abs(fp - np.median(fp)) > 3), 2),
        tau_eq22_mean_ns=round(tau.mean() * 1e9, 1), tau_eq22_std_ns=round(tau.std() * 1e9, 1), tau_eq22_max_abs_ns=round(np.abs(tau).max() * 1e9, 0),
        psi_resultant_length=round(R_psi, 4), psi_lag1_wrapped_diff_std_rad=round(np.angle(np.exp(1j * np.diff(psi))).std(), 3),
        zero_rows="|".join(map(str, np.where(pw_tone < 1e-9 * pw_tone.max())[0])), notes=c["notes"]))
    prof[f"{c['cls']}_{c['capture']}_{c['date']}"] = 10 * np.log10(pw_tone + 1e-12)
    print(c["cls"], c["date"], c["capture"], n, rows[-1]["duration_s"], "s", flush=True)
df = pd.DataFrame(rows); df.to_csv(os.path.join(a.out, "dataset_inventory.csv"), index=False)
pd.DataFrame(prof).assign(row=np.arange(256)).to_csv(os.path.join(a.out, "tone_power_profile.csv"), index=False)
print(df.drop(columns=["notes", "code"]).to_string())
