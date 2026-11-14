"""Step 2: real-data KPI experiment on the Brescia 20 MHz CSI (protocol: docs/protocol_m4.md).
Per capture: skip first 200 frames, cut non-overlapping windows of P=300 frames (discard windows whose pcap timestamps contain a gap > 0.5 s),
apply raw_nogain / A0 (Eq.6) / A1 (Eq.6+Eq.21) / A2 (Eq.6+Eq.22) / A3 (Eq.6+Algorithm 4+Eq.25) and compute KPI1-KPI4 + diagnostics.
Outputs: windows.csv (one row per capture/window/method), windows_diag.csv (per capture/window diagnostics), run_info.json.
Usage: python 02_real_kpi.py --root /path/to/20MHz --out ../Results/real_data/run1 [--maxwin 0] [--procs 4]"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import csi_methods as M, load_brescia as L

METHODS = ["raw_nogain", "A0", "A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"]
KIDX, TS, P, SKIP, GAP = M.HE20_KIDX, M.TS_HE, 300, 200, 0.5

def vdb(X): return M.dynamic_power_fraction_db(X)

def one_capture(c):
    X = L.load_active(c["mat"])[SKIP:]; ts = L.pcap_timestamps(c["pcap"])[SKIP:]
    dt = np.diff(ts); nwin = X.shape[0] // P
    rows, diag = [], []
    for w in range(nwin):
        sl = slice(w * P, (w + 1) * P)
        if np.any(dt[sl.start:sl.stop - 1] > GAP): continue
        Xw = X[sl]; h = P // 2; base = dict(cls=c["cls"], people=c["people"], date=c["date"], capture=c["capture"], window=w, t_start_s=round(ts[sl.start] - ts[0], 1))
        G = M.gain_eq6(Xw)
        # --- per-method KPIs (timing covers Eq.6 + estimator + apply, as in M.correct) ---
        Ys, TIMES = {}, {}
        for m in METHODS:
            t0 = time.perf_counter(); Y = M.correct_m4(m, Xw, KIDX, TS); ms = (time.perf_counter() - t0) * 1e3; Ys[m] = Y; TIMES[m] = ms
        mask = M.static_tone_mask(M.gain_eq6(Ys["A3_algo4"]))      # same tone set for every method (chosen from the Algorithm 4 output)
        for m in METHODS:
            Y = Ys[m]; cm, cl = M.complex_consistency(Y)
            rows.append(dict(**base, method=m, V_dB=vdb(Y), phase_circstd_deg=M.phase_circ_std_deg_masked(M.gain_eq6(Y), mask), C_mean=cm, C_lag1=cl, time_ms=TIMES[m],
                             V2_insample_dB=vdb(Y[h:]), V2_oos_dB=np.nan))
        # --- Algorithm 4 out-of-sample: reference from the first half only, scored on the second half ---
        tc, pc = M.eq22(G[:h], KIDX, TS); b_first = M.apply_phase(G[:h], tc, pc, KIDX, TS).mean(axis=0)
        rows[-1]["V2_oos_dB"] = vdb(M.correct_m4("A3_algo4", Xw[h:], KIDX, TS, b_ref=b_first))
        # --- Eq.(21)/(22)/Algorithm 4 validity diagnostics on the same window ---
        t1, p1 = M.eq21(G, KIDX, TS); t1w, _ = M.eq21_wide(G, KIDX, TS); t2, p2 = M.eq22(G, KIDX, TS); t3, p3 = M.algo4_eq25(G, KIDX, TS)
        tc2, pc2 = M.eq22(G, KIDX, TS); b = M.apply_phase(G, tc2, pc2, KIDX, TS).mean(axis=0)
        kb = np.abs(b) ** 2 > 0.1
        d21 = (t1 - np.median(t1)) - (t3 - np.median(t3)); d21w = (t1w - np.median(t1w)) - (t3 - np.median(t3)); d22 = (t2 - np.median(t2)) - (t3 - np.median(t3))
        fp = 10 * np.log10(np.mean(np.abs(Xw) ** 2, axis=1))
        diag.append(dict(**base, frame_power_std_dB=fp.std(), frame_power_range_dB=fp.max() - fp.min(),
            tau21_med_ns=np.median(t1) * 1e9, tau22_med_ns=np.median(t2) * 1e9, tau3_med_ns=np.median(t3) * 1e9,
            tau21_std_ns=t1.std() * 1e9, tau22_std_ns=t2.std() * 1e9, tau3_std_ns=t3.std() * 1e9, tau22_absmax_ns=np.abs(t2).max() * 1e9,
            tau22_over_halfTs=float(np.mean(np.abs(t2) > TS / 2)), kbar_frac=kb.mean(), b_flatness_std_dB=np.std(10 * np.log10(np.abs(b) ** 2 + 1e-12)),
            d21_vs_A3_med_abs_ns=np.median(np.abs(d21)) * 1e9, d21_vs_A3_frac_gt25ns=float(np.mean(np.abs(d21) > 25e-9)),
            d21w_vs_A3_med_abs_ns=np.median(np.abs(d21w)) * 1e9, d21w_vs_A3_frac_gt25ns=float(np.mean(np.abs(d21w) > 25e-9)), d22_vs_A3_med_abs_ns=np.median(np.abs(d22)) * 1e9, d22_vs_A3_frac_gt25ns=float(np.mean(np.abs(d22) > 25e-9))))
    print(f"{c['cls']:8s} {c['date']} {c['capture']}: windows {len(diag)}/{nwin}", flush=True)
    return rows, diag

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--maxcap", type=int, default=0, help="debug: only first N captures"); ap.add_argument("--procs", type=int, default=4); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); caps = L.list_captures(a.root)
    if a.maxcap: caps = caps[:a.maxcap]
    t0 = time.time()
    with Pool(a.procs) as pool: res = pool.map(one_capture, caps)
    rows = [r for rr, _ in res for r in rr]; diag = [d for _, dd in res for d in dd]
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "windows.csv"), index=False); pd.DataFrame(diag).to_csv(os.path.join(a.out, "windows_diag.csv"), index=False)
    json.dump(dict(script="02_real_kpi.py", root=a.root, P=P, skip=SKIP, gap_s=GAP, tones=len(KIDX), Ts=TS, n_captures=len(caps), n_rows=len(rows), seconds=round(time.time() - t0, 1),
                   note="independent unit = capture; windows within a capture are correlated"), open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    print("done", round(time.time() - t0, 1), "s")
