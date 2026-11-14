"""Step 3: aggregate windows.csv -> capture-level and class-level tables, protocol pass/fail rules, occupancy-preservation metrics,
Eq.(21)/(22) validity diagnostics.  Independent unit = capture (mean over its windows). Usage: python 03_analyze.py --run ../Results/real_data/run1"""
import argparse, os, numpy as np, pandas as pd
from itertools import product

ap = argparse.ArgumentParser(); ap.add_argument("--run", required=True); a = ap.parse_args()
T = os.path.join(a.run, "tables"); os.makedirs(T, exist_ok=True)
W = pd.read_csv(os.path.join(a.run, "windows.csv")); D = pd.read_csv(os.path.join(a.run, "windows_diag.csv"))
inv = pd.read_csv(os.path.join(a.run, "../../../Data/dataset_inventory.csv"))
full = set(inv.loc[inv.full_length == 1, ["cls", "date", "capture"]].apply(tuple, axis=1))
W["full"] = [(c, d, k) in full for c, d, k in zip(W.cls, W.date, W.capture)]; D["full"] = [(c, d, k) in full for c, d, k in zip(D.cls, D.date, D.capture)]
M = ["raw_nogain", "A0", "A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"]; CL = ["E_day", "E_night", "S1", "M4"]
out = []; P = lambda *s: out.append(" ".join(str(x) for x in s))

# ---- capture-level means (mean over windows) ----
cap = W.groupby(["cls", "date", "capture", "full", "method"]).agg(n_win=("window", "count"), V_dB=("V_dB", "mean"), V_std_over_windows=("V_dB", "std"), phase_circstd_deg=("phase_circstd_deg", "mean"),
      C_mean=("C_mean", "mean"), C_lag1=("C_lag1", "mean"), time_ms=("time_ms", "mean"), V2_in=("V2_insample_dB", "mean"), V2_oos=("V2_oos_dB", "mean")).reset_index()
cap.to_csv(os.path.join(T, "capture_level_kpis.csv"), index=False)

# ---- class-level (primary = full captures only) ----
def classtab(df, cols):
    g = df.groupby(["cls", "method"])
    res = []
    for (c, m), x in g:
        r = dict(cls=c, method=m, n_captures=len(x))
        for col in cols: r[col + "_mean"] = x[col].mean(); r[col + "_min"] = x[col].min(); r[col + "_max"] = x[col].max()
        res.append(r)
    r = pd.DataFrame(res); r["method"] = pd.Categorical(r.method, M); r["cls"] = pd.Categorical(r.cls, CL)
    return r.sort_values(["cls", "method"]).reset_index(drop=True)
prim = cap[cap.full]; sec = cap[~cap.full]
ct = classtab(prim, ["V_dB", "phase_circstd_deg", "C_mean", "C_lag1", "time_ms"]); ct.to_csv(os.path.join(T, "class_level_kpis_primary.csv"), index=False)
classtab(sec, ["V_dB", "phase_circstd_deg", "C_mean", "C_lag1"]).to_csv(os.path.join(T, "class_level_kpis_secondary_partial_captures.csv"), index=False)
P("=== KPI1  V (dB, lower better): mean over captures [min, max], PRIMARY (full-length captures) ===")
for c in CL:
    P(f"{c:8s} n={int(ct[ct.cls==c].n_captures.iloc[0])}  " + "  ".join(f"{m}: {r.V_dB_mean:7.2f} [{r.V_dB_min:7.2f},{r.V_dB_max:7.2f}]" for m, r in ct[ct.cls == c].set_index('method').iterrows()))
P("\n=== KPI2 median phase circ-std (deg), KPI4 C_mean / C_lag1, KPI3 ms/window (primary) ===")
for c in CL:
    for m, r in ct[ct.cls == c].set_index("method").iterrows():
        P(f"{c:8s} {m:11s} phase {r.phase_circstd_deg_mean:7.2f}  C_mean {r.C_mean_mean:6.4f}  C_lag1 {r.C_lag1_mean:6.4f}  time {r.time_ms_mean:6.2f} ms")

# ---- Rule (a) ----
P("\n=== Rule (a): V(A3) < V(A0), V(A1), V(A2) in EVERY full capture of E_day, E_night, S1 ===")
ra = []
for (c, d, k), x in prim[prim.cls.isin(["E_day", "E_night", "S1"])].groupby(["cls", "date", "capture"]):
    v = x.set_index("method").V_dB; ok = all(v["A3_algo4"] < v[m] for m in ["A0", "A1_eq21", "A2_eq22"]); ok_w = ok and v["A3_algo4"] < v["A1w_eq21_wide"]
    ra.append(dict(cls=c, date=d, capture=k, A0=v["A0"], A1=v["A1_eq21"], A1w=v["A1w_eq21_wide"], A2=v["A2_eq22"], A3=v["A3_algo4"], A3_minus_A2=v["A3_algo4"] - v["A2_eq22"], A3_minus_A1=v["A3_algo4"] - v["A1_eq21"], A3_minus_A1w=v["A3_algo4"] - v["A1w_eq21_wide"], rule_a_holds_vs_A1_as_protocol=ok, rule_a_holds_vs_A1w_strict=ok_w))
ra = pd.DataFrame(ra); ra.to_csv(os.path.join(T, "rule_a_per_capture.csv"), index=False); P(ra.round(2).to_string(index=False))
P("RULE (a) overall vs A1 as frozen (M3 code, n in [-3,3]):", "PASS" if ra.rule_a_holds_vs_A1_as_protocol.all() else f"FAIL ({(~ra.rule_a_holds_vs_A1_as_protocol).sum()} of {len(ra)} captures violate)")
P("RULE (a) overall vs A1w (widened n in [-8,8], stricter, added post hoc):", "PASS" if ra.rule_a_holds_vs_A1w_strict.all() else f"FAIL ({(~ra.rule_a_holds_vs_A1w_strict).sum()} of {len(ra)} captures violate)")
# window-level descriptive: how often does A3 beat A2 / A1 within a window
wp = W[W.full].pivot_table(index=["cls", "date", "capture", "window"], columns="method", values="V_dB").reset_index()
wp["A3_lt_A2"] = wp.A3_algo4 < wp.A2_eq22; wp["A3_lt_A1"] = wp.A3_algo4 < wp.A1_eq21; wp["A3_lt_A1w"] = wp.A3_algo4 < wp.A1w_eq21_wide
P("\nWindow-level share where A3 has lower V than A2 / A1 (descriptive, windows correlated):")
P(wp.groupby("cls")[["A3_lt_A2", "A3_lt_A1", "A3_lt_A1w"]].mean().round(3).to_string())

# ---- Rule (b) ----
P("\n=== Rule (b): |V2_oos - V2_in| for Algorithm 4 on matched frames (second half of window), class mean over primary captures ===")
rb = []
for c in CL:
    x = prim[(prim.cls == c) & (prim.method == "A3_algo4")]
    rb.append(dict(cls=c, V2_in=x.V2_in.mean(), V2_oos=x.V2_oos.mean(), diff_dB=(x.V2_oos - x.V2_in).mean(), max_abs_diff_over_captures=(x.V2_oos - x.V2_in).abs().max(), rule_b_holds=abs((x.V2_oos - x.V2_in).mean()) < 1))
rb = pd.DataFrame(rb); rb.to_csv(os.path.join(T, "rule_b_in_vs_oos.csv"), index=False); P(rb.round(3).to_string(index=False))
P("RULE (b) overall:", "PASS" if rb.rule_b_holds.all() else "FAIL (see per-class)")

# ---- Rule (c') and occupancy preservation ----
P("\n=== Rule (c') / occupancy preservation: same-day pairs. separation = min(V occupied captures) - max(V empty captures), dB (>0 = every occupied capture above every empty capture) ===")
pairs = [("S1", "E_day"), ("M4", "E_night")]; rc = []
for m in M:
    for occ, emp in pairs:
        vo = prim[(prim.cls == occ) & (prim.method == m)].V_dB; ve = prim[(prim.cls == emp) & (prim.method == m)].V_dB
        # window-level AUC (descriptive)
        wo = W[(W.full) & (W.cls == occ) & (W.method == m)].V_dB.values; we = W[(W.full) & (W.cls == emp) & (W.method == m)].V_dB.values
        ranks = pd.Series(np.r_[wo, we]).rank().values; auc = (ranks[:len(wo)].sum() - len(wo) * (len(wo) + 1) / 2) / (len(wo) * len(we))
        rc.append(dict(method=m, pair=f"{occ} vs {emp}", V_occ_mean=vo.mean(), V_emp_mean=ve.mean(), delta_V_dB=vo.mean() - ve.mean(), separation_dB=vo.min() - ve.max(), window_AUC_occ_above_emp=auc, n_occ=len(vo), n_emp=len(ve)))
rc = pd.DataFrame(rc); rc.to_csv(os.path.join(T, "rule_c_occupancy_preservation.csv"), index=False); P(rc.round(3).to_string(index=False))
a3 = rc[rc.method == "A3_algo4"]; P("RULE (c') for A3:", "; ".join(f"{r.pair}: separation {r.separation_dB:.2f} dB -> {'PASS' if r.separation_dB > 0 else 'FAIL'}" for r in a3.itertuples()))

# ---- power-domain activity (what Eq.6 removes) ----
P("\n=== Frame-power std before Eq.(6) (dB), per class: mean over primary captures [min,max] and share of windows with a >3 dB frame-power range ===")
dg = D[D.full].groupby(["cls", "capture"]).agg(fp_std=("frame_power_std_dB", "mean"), fp_range=("frame_power_range_dB", "mean"), gt3=("frame_power_range_dB", lambda x: (x > 3).mean())).reset_index()
dg.to_csv(os.path.join(T, "frame_power_by_capture.csv"), index=False)
for c in CL:
    x = dg[dg.cls == c]; P(f"{c:8s} std {x.fp_std.mean():.3f} [{x.fp_std.min():.3f},{x.fp_std.max():.3f}] dB ; windows with range>3 dB: {x.gt3.mean()*100:.1f}%")

# ---- validity of Eq.(21)/(22)/Alg.4 ----
P("\n=== Validity diagnostics (primary), mean over captures of window means ===")
dv = D[D.full].groupby(["cls", "capture"]).mean(numeric_only=True).reset_index().groupby("cls").mean(numeric_only=True)
cols = ["tau22_med_ns", "tau22_std_ns", "tau22_absmax_ns", "tau22_over_halfTs", "kbar_frac", "b_flatness_std_dB", "d21_vs_A3_med_abs_ns", "d21_vs_A3_frac_gt25ns", "d21w_vs_A3_med_abs_ns", "d21w_vs_A3_frac_gt25ns", "d22_vs_A3_med_abs_ns", "d22_vs_A3_frac_gt25ns", "tau21_std_ns", "tau3_std_ns"]
dv[cols].to_csv(os.path.join(T, "validity_diagnostics_by_class.csv")); P(dv.loc[CL, cols].round(3).T.to_string())
# static-fraction proxy
gam = prim[prim.method == "A3_algo4"].groupby("cls").V_dB.mean(); P("\nStatic-fraction proxy gamma_hat = 1 - 10^(V_A3/10) per class:", {c: round(1 - 10 ** (gam[c] / 10), 5) for c in CL})
open(os.path.join(a.run, "summary.txt"), "w").write("\n".join(out)); print("\n".join(out))
