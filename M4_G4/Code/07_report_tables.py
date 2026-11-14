"""Step 7: turn the result CSVs into Report/report_data.json (strings ready for the report tables) so no number is typed by hand."""
import json, numpy as np, pandas as pd
R = "../Results/real_data/run1/tables/"; f2 = lambda x, d=2: f"{x:.{d}f}"
CL = ["E_day", "E_night", "S1", "M4"]; MM = ["A0", "A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"]
ct = pd.read_csv(R + "class_level_kpis_primary.csv"); inv = pd.read_csv("../Data/dataset_inventory.csv"); out = {}
g = lambda c, m: ct[(ct.cls == c) & (ct.method == m)].iloc[0]
out["kpi1"] = [[c, str(int(g(c, "A0").n_captures))] + [("%s" % f2(g(c, m).V_dB_mean)) if m == "A0" else f"{f2(g(c, m).V_dB_mean)} [{f2(g(c, m).V_dB_min)}, {f2(g(c, m).V_dB_max)}]" for m in MM] for c in CL]
out["kpi24"] = [[c, m.split("_")[0], f2(g(c, m).phase_circstd_deg_mean, 2), f2(g(c, m).C_mean_mean, 4), f2(g(c, m).C_lag1_mean, 4)] for c in CL for m in MM]
out["time"] = [[m.split("_")[0], f2(ct[ct.method == m].time_ms_mean.min(), 1), f2(ct[ct.method == m].time_ms_mean.max(), 1)] for m in MM]
ra = pd.read_csv(R + "rule_a_per_capture.csv"); out["rule_a"] = [[r.cls, r.capture.replace("capture", "c"), f2(r.A1), f2(r.A1w), f2(r.A2), f2(r.A3), f2(r.A3_minus_A1w), f2(r.A3_minus_A2), "yes" if r.rule_a_holds_vs_A1w_strict else "NO"] for r in ra.itertuples()]
rb = pd.read_csv(R + "rule_b_in_vs_oos.csv"); out["rule_b"] = [[r.cls, f2(r.V2_in), f2(r.V2_oos), f2(r.diff_dB, 3), f2(r.max_abs_diff_over_captures, 3)] for r in rb.itertuples()]
rc = pd.read_csv(R + "rule_c_occupancy_preservation.csv"); rc = rc[rc.method.isin(["A1_eq21", "A1w_eq21_wide", "A2_eq22", "A3_algo4"])]
out["occ"] = [[r.method.split("_")[0], r.pair, f2(r.V_occ_mean), f2(r.V_emp_mean), f2(r.delta_V_dB), f2(r.separation_dB), f2(r.window_AUC_occ_above_emp, 3)] for r in rc.itertuples()]
gs = pd.read_csv("../Data/gain_structure.csv"); gg = gs.groupby("cls").agg(j=("jumps_per_1000_frames", "mean"), jmin=("jumps_per_1000_frames", "min"), jmax=("jumps_per_1000_frames", "max"), big=("median_big_jump_dB", "median"), sa=("pow_std_all_dB", "mean"), sm=("pow_std_within_mid_state_dB", "mean"), smin=("pow_std_within_mid_state_dB", "min"), smax=("pow_std_within_mid_state_dB", "max")).reindex(CL)
out["gain"] = [[c, f2(r.j, 1), f"{f2(r.jmin,1)}-{f2(r.jmax,1)}", f2(r.big, 2), f2(r.sa, 2), f"{f2(r.sm,2)} ({f2(r.smin,2)}-{f2(r.smax,2)})"] for c, r in gg.iterrows()]
fl = pd.read_csv(R + "eq21_failures.csv"); ff = fl[fl.frames > 5000].groupby("cls").agg(pm=("pct", "mean"), pmin=("pct", "min"), pmax=("pct", "max"), n=("frames_eq21_vs_A3_tau_gt25ns", "sum")).reindex(CL)
out["eq21fail"] = [[c, f2(r.pm, 3), f"{f2(r.pmin,3)}-{f2(r.pmax,3)}", str(int(r.n))] for c, r in ff.iterrows()]
vd = pd.read_csv(R + "validity_diagnostics_by_class.csv", index_col=0).reindex(CL)
out["valid"] = [[c, f2(r.tau22_med_ns, 0), f2(r.tau22_std_ns, 0), f2(r.tau22_absmax_ns, 0), f2(100 * r.tau22_over_halfTs, 1), f2(r.kbar_frac, 3), f2(r.b_flatness_std_dB, 1)] for c, r in vd.iterrows()]
out["inv"] = [[r.cls, r.date, r.capture.replace("capture", "c"), str(r.frames_mat), f2(r.duration_s / 60, 1), f2(r.dt_median_ms, 1), str(r.n_gaps_gt_0p5s) + (f" (max {r.max_gap_s:.0f} s)" if r.n_gaps_gt_0p5s else ""), str(r.usable_windows_P300), "primary" if r.full_length else "partial"] for r in inv.itertuples()]
S = pd.read_csv("../Results/synthetic/SYNTHETIC_known_answer_V.csv"); piv = S.pivot_table(index=["channel", "gamma"], columns="method", values="V_dB_mean").reset_index(); fl_ = S[S.method == "A3_algo4"].set_index(["channel", "gamma"]).known_floor_dB
out["synth"] = [[r.channel, str(r.gamma), f2(fl_[(r.channel, r.gamma)]), f2(r.A1_eq21), f2(r.A1w_eq21_wide), f2(r.A2_eq22), f2(r.A3_algo4)] for r in piv.itertuples()]
sec = pd.read_csv(R + "class_level_kpis_secondary_partial_captures.csv"); out["sec"] = [[r.cls, m.split("_")[0], f2(r.V_dB_mean)] for r in sec.itertuples() for m in [r.method] if m in MM]
json.dump(out, open("../Report/report_data.json", "w"), indent=1); print({k: len(v) for k, v in out.items()}); print(out["gain"]); print(out["eq21fail"]); print(out["time"])
