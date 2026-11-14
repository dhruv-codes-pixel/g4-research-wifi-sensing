"""Step 4 (diagnostics, REAL data): (1) gain-state structure of every capture (what Eq.(6) removes), (2) Eq.(21) failure analysis
(2*pi integer ambiguity across the 3-tone DC gap).  Outputs: Data/gain_structure.csv, Results/real_data/run1/tables/eq21_failures.csv,
eq21_ambiguity_examples.csv.  Usage: python 04_diagnostics.py --root /path/to/20MHz"""
import argparse, os, numpy as np, pandas as pd
import csi_methods as M, load_brescia as L
K, TS = M.HE20_KIDX, M.TS_HE
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--res", default="../Results/real_data/run1/tables"); ap.add_argument("--data", default="../Data"); a = ap.parse_args()
os.makedirs(a.res, exist_ok=True)

# ---------- (1) gain-state structure ----------
rows = []
for c in L.list_captures(a.root):
    X = L.load_active(c["mat"])[200:]; fp = 10 * np.log10(np.mean(np.abs(X) ** 2, axis=1)); d = np.diff(fp); med = np.median(fp)
    big = np.abs(d) > 3; mid = np.abs(fp - med) < 3
    # flat-ratio test of the big jumps: a scalar gain step changes every tone by the same dB
    hi = fp > med + 3
    flat = np.nan
    if hi.sum() > 50 and (mid & (np.abs(fp - med) < 0.7)).sum() > 50:
        r = 20 * np.log10(np.mean(np.abs(X[hi]), axis=0) / np.mean(np.abs(X[mid & (np.abs(fp - med) < 0.7)]), axis=0)); flat = r.std()
    rows.append(dict(cls=c["cls"], date=c["date"], capture=c["capture"], frames=len(fp), jumps_gt3dB=int(big.sum()), jumps_per_1000_frames=round(1000 * big.mean(), 2),
        median_big_jump_dB=round(float(np.median(np.abs(d[big]))), 2) if big.any() else np.nan, share_frames_gt3dB_above_median_pct=round(100 * hi.mean(), 2),
        pow_std_all_dB=round(fp.std(), 3), pow_std_within_mid_state_dB=round(fp[mid].std(), 3), p_change_gt0p15dB=round(float(np.mean(np.abs(d) > 0.15)), 3),
        tone_ratio_std_dB_hi_vs_mid=round(flat, 3) if flat == flat else np.nan))
g = pd.DataFrame(rows); g.to_csv(os.path.join(a.data, "gain_structure.csv"), index=False); print(g.to_string(index=False))

# ---------- (2) Eq.(21) failures vs Algorithm 4 (diagnostic only; A3 used as reference for the comparison) ----------
def eq21_forced(G, n):
    """Eq.(21) with the 2*pi offset of the second tone block forced to n (copy of the M3 logic for a single n)."""
    f = K / TS; ph = np.angle(G); U = np.empty_like(ph)
    for b in M.blocks(K): U[:, b] = np.unwrap(ph[:, b], axis=1)
    bl = M.blocks(K); Y = U.copy(); Y[:, bl[1]] += 2 * np.pi * n
    A = np.stack([-2 * np.pi * f, -np.ones_like(f)], axis=1); th = Y @ np.linalg.pinv(A).T
    return th[:, 0], th[:, 1], np.sum((Y - th @ A.T) ** 2, axis=1)
fr, ex = [], []
for c in L.list_captures(a.root):
    X = L.load_active(c["mat"])[200:]; nb = 0; nt = 0; tot_bad_d = []
    for w in range(X.shape[0] // 300):
        Xw = X[w * 300:(w + 1) * 300]; G = M.gain_eq6(Xw); t1, p1 = M.eq21(G, K, TS); t3, p3 = M.algo4_eq25(G, K, TS)
        d = (t1 - np.median(t1)) - (t3 - np.median(t3)); bad = np.abs(d) > 25e-9; nb += bad.sum(); nt += len(d); tot_bad_d += list(d[bad] * 1e9)
        if bad.any() and len(ex) < 6 and c["cls"] in ("E_day", "S1"):
            i = int(np.where(bad)[0][0]); res = []
            for n in range(-3, 4):
                t, p, sse = eq21_forced(G[i:i + 1], n); res.append((n, t[0] * 1e9, sse[0]))
            n_best = min(res, key=lambda r: r[2])[0]
            ex.append(dict(cls=c["cls"], capture=c["capture"], window=w, frame=i, tau_A3_minus_window_median_ns=(t3[i] - np.median(t3)) * 1e9, tau_eq21_minus_window_median_ns=(t1[i] - np.median(t1)) * 1e9,
                           n_chosen_by_eq21=n_best, sse_by_n={n: round(s, 2) for n, _, s in res}, tau_by_n_ns={n: round(t, 1) for n, t, _ in res}))
    fr.append(dict(cls=c["cls"], date=c["date"], capture=c["capture"], frames=nt, frames_eq21_vs_A3_tau_gt25ns=int(nb), pct=round(100 * nb / nt, 3),
                   median_abs_d_bad_ns=round(float(np.median(np.abs(tot_bad_d))), 1) if tot_bad_d else np.nan))
    print(fr[-1], flush=True)
pd.DataFrame(fr).to_csv(os.path.join(a.res, "eq21_failures.csv"), index=False); pd.DataFrame(ex).to_csv(os.path.join(a.res, "eq21_ambiguity_examples.csv"), index=False)
print(pd.DataFrame(ex).to_string(index=False))
