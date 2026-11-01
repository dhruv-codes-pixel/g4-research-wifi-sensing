"""
make_gain_figures.py -- builds the M3 gain figures/tables from the saved raw data (no new simulation).
Inputs : ../Results/gain/{E2,E5,FIG6}.npz, fig6_paper_digitised.csv, sanity_test.csv
Outputs: ../Results/gain/figure6_comparison.png, lambda_analysis.png, before_after.csv, fig6_table.csv
Run    : python make_gain_figures.py
"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csi_sim as c
import gain_fix as gf

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Results", "gain")
L = lambda f: np.load(os.path.join(R, f), allow_pickle=False)


def med_lin(chi, paper=True):
    v = chi ** 2 / (1 - chi ** 2) if paper else chi / (1 - chi)
    return np.nanmedian(v, axis=-1)


def db(x):
    x = np.asarray(x, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(x > 0, 10 * np.log10(np.where(x > 0, x, 1)), np.nan)


# ------------------------------------------------------------------ Figure 6 comparison
def figure6():
    z = L("FIG6.npz"); x = z["x"]; meths = list(z["methods"])
    paper = {}
    for r in csv.DictReader(open(os.path.join(R, "fig6_paper_digitised.csv"))):
        paper.setdefault((r["dynamic_type"], r["series"]), []).append((float(r["gamma"]), float(r["snr_metric_chi2_over_1_minus_chi2"])))
    pairs = [("Baseline (6)", "Eq6 (baseline)", "tab:red"), ("Baseline (8)", "Eq8 (baseline)", "tab:green"),
             ("Algorithm 1", "Algo1", "tab:blue"), ("Algorithm 3", "Algo3 fixed", "black")]
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 9.5)); rows = []
    for col, (dyn, tag) in enumerate([("i", "paper_i"), ("ii", "paper_ii")]):
        chi = z[f"{tag}__chi"]
        ours = {m: med_lin(chi[i]) for i, m in enumerate(meths)}
        a = ax[0, col]
        for pn, on, cl in pairs:
            if (dyn, pn) in paper:
                g, v = np.array(paper[(dyn, pn)]).T; a.plot(g, v, ":", marker="s", mfc="none", color=cl, label=f"paper {pn}")
            a.plot(x, ours[on], "-", marker="o", ms=4, color=cl, label=f"ours {on}")
        a.plot(x, ours["Algo3 M2-original"], "-", marker="^", ms=4, color="tab:gray", label="ours Algo3 M2-original")
        a.plot(x, ours["Ideal gain (oracle)"], "--", color="tab:purple", lw=2, label="ours ideal gain (ceiling)")
        a.set_yscale("log"); a.set_title(f"Fig. 6({'ab'[col]}) dynamic type ({dyn}), K=256, 500 trials", fontsize=10)
        a.set_ylabel("median  chi^2/(1-chi^2)"); a.grid(alpha=.3, which="both"); a.legend(fontsize=6.5, ncol=2)
        b = ax[1, col]                                 # relative to baseline (6): removes the absolute-scale difference
        e6p = dict(paper.get((dyn, "Baseline (6)"), []))
        for pn, on, cl in pairs[1:]:
            if (dyn, pn) in paper:
                pts = [(g, v / e6p[g]) for g, v in paper[(dyn, pn)] if g in e6p]
                if pts: g, v = np.array(pts).T; b.plot(g, db(v), ":", marker="s", mfc="none", color=cl, label=f"paper {pn} / (6)")
            b.plot(x, db(ours[on] / ours["Eq6 (baseline)"]), "-", marker="o", ms=4, color=cl, label=f"ours {on} / Eq6")
        b.plot(x, db(ours["Ideal gain (oracle)"] / ours["Eq6 (baseline)"]), "--", color="tab:purple", lw=2, label="ours ideal / Eq6")
        b.axhline(0, color="k", lw=.8); b.set_ylabel("ratio to Baseline (6)  (dB)"); b.set_xlabel("gamma")
        b.grid(alpha=.3); b.legend(fontsize=6.5, ncol=2)
        for j, g in enumerate(x):
            row = dict(dynamic_type=dyn, gamma=g)
            for pn, on, _ in pairs:
                pv = dict(paper.get((dyn, pn), [])).get(round(g, 2))
                row[f"paper_{pn}"] = None if pv is None else round(pv, 2); row[f"ours_{on}"] = round(float(ours[on][j]), 2)
            row["ours_Algo3 M2-original"] = round(float(ours["Algo3 M2-original"][j]), 2)
            row["ours_Ideal gain (oracle)"] = round(float(ours["Ideal gain (oracle)"][j]), 2)
            rows.append(row)
    plt.tight_layout(); plt.savefig(os.path.join(R, "figure6_comparison.png"), dpi=150); plt.close()
    with open(os.path.join(R, "fig6_table.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return rows


# ------------------------------------------------------------------ lambda analysis
def lambda_analysis():
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 9))
    s = c.simulate(K=256, gamma=0.9, rng=np.random.default_rng([42, 5, 2, 0]))   # = E5, M2 simulator, K=256, trial 0 (not selected)
    _, io = gf.gain_algo3(s["ht"], "original"); _, ifx = gf.gain_algo3(s["ht"], "fixed")
    To, Tf = io["table"], ifx["table"]
    a = ax[0, 0]
    a.plot(Tf[:, 0], Tf[:, 1], "-", color="tab:gray", lw=1, label="mean(Gh^2), fine grid")
    a.plot(Tf[:, 0], Tf[:, 2], "--", color="tab:red", label="paper test threshold lambda^2/24")
    a.plot(To[:, 0], To[:, 1], "o", color="tab:green", label="M2 grid points (0.05 lambda_max steps)")
    rej = To[:, 4] == 0; a.plot(To[rej, 0], To[rej, 1], "x", color="k", ms=8, label="rejected by lambda^2/24")
    a.axvline(0.5, color="tab:purple", ls=":", label="true AGC step 0.5 dB")
    a.axvline(io["lam_hat"], color="tab:green", ls="-.", label=f"M2 choice {io['lam_hat']:.2f} dB")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("lambda (dB)"); a.set_title("Why M2 picks ~3 dB (E5, M2 simulator, K=256, trial 0)", fontsize=10)
    a.legend(fontsize=7); a.grid(alpha=.3, which="both")
    a = ax[0, 1]
    P = 300; thr = np.sqrt(np.log(1e3) / P)
    a.plot(Tf[:, 0], Tf[:, 3], "-", color="tab:gray", label="m = |mean exp(j2pi Gh/lambda)|")
    a.axhline(thr, color="tab:purple", ls="--", label="FIX-2 Rayleigh threshold (alpha=1e-3)")
    obj = np.where(np.isfinite(Tf[:, 6]), Tf[:, 6], np.nan)
    a2 = a.twinx(); a2.plot(Tf[:, 0], obj, ".", ms=3, color="tab:orange", label="Obj(lambda), accepted"); a2.set_yscale("log"); a2.set_ylabel("Obj")
    a.axvline(0.5, color="tab:purple", ls=":"); a.axvline(ifx["lam_hat"], color="tab:orange", ls="-.", label=f"fixed choice {ifx['lam_hat']:.2f} dB")
    a.set_xscale("log"); a.set_xlabel("lambda (dB)"); a.set_ylabel("m"); a.legend(fontsize=7, loc="upper left"); a2.legend(fontsize=7, loc="upper right")
    a.set_title("Fixed: accept by Eq.18 statistic, fine grid (same trial)", fontsize=10)
    # distributions of lambda_hat (E5, both simulators)
    a = ax[1, 0]; z = L("E5.npz"); xs = z["x"]
    pos = np.arange(len(xs)); w = 0.18
    for k, (tag, row, cl, lab) in enumerate([("m2", 0, "tab:green", "M2 sim, original"), ("m2", 1, "tab:purple", "M2 sim, fixed"),
                                             ("paper", 0, "tab:olive", "paper g1, original"), ("paper", 1, "black", "paper g1, fixed")]):
        lam = z[f"{tag}__lam"][row]
        a.boxplot([lam[j][np.isfinite(lam[j])] for j in range(len(xs))], positions=pos + (k - 1.5) * w, widths=w * .9,
                  showfliers=False, patch_artist=True, boxprops=dict(facecolor=cl, alpha=.35), medianprops=dict(color=cl))
        a.plot([], [], color=cl, lw=6, alpha=.4, label=lab)
    a.axhline(0.5, color="tab:purple", ls=":"); a.set_xticks(pos); a.set_xticklabels([f"K={int(v)}" for v in xs])
    a.set_yscale("log"); a.set_ylabel("selected lambda (dB)"); a.set_title("E5 (gamma=0.9): selected lambda, 500 trials", fontsize=10)
    a.legend(fontsize=7); a.grid(alpha=.3, which="both")
    # lambda vs noise from the unit tests
    a = ax[1, 1]; rows = list(csv.DictReader(open(os.path.join(R, "sanity_test.csv"))))
    for drift, ls in (("False", "-"), ("True", "--")):
        for var, cl in (("original", "tab:green"), ("fixed", "tab:purple"), ("fixed+bonf", "tab:brown")):
            rr = [r for r in rows if r["large_scale_drift"] == drift and r["variant"] == var]
            nx = [float(r["noise_dB"]) for r in rr]; md = [float(r["est_lambda_median_dB"]) for r in rr]
            lo = [float(r["est_lambda_q25"]) for r in rr]; hi = [float(r["est_lambda_q75"]) for r in rr]
            a.errorbar(np.array(nx) + (0.004 if drift == "True" else 0), md, yerr=[np.subtract(md, lo), np.subtract(hi, md)],
                       fmt=ls + "o", ms=4, color=cl, capsize=2, label=f"{var}{' + slow drift' if drift == 'True' else ''}")
    a.axhline(0.5, color="tab:purple", ls=":"); a.set_yscale("log"); a.set_xlabel("noise = std of true power fluctuation (dB)")
    a.set_ylabel("selected lambda (dB), median & IQR"); a.set_title("Unit tests: lambda vs noise (100 trials each)", fontsize=10)
    a.legend(fontsize=6.5); a.grid(alpha=.3, which="both")
    plt.tight_layout(); plt.savefig(os.path.join(R, "lambda_analysis.png"), dpi=150); plt.close()
    return io["lam_hat"], ifx["lam_hat"]


# ------------------------------------------------------------------ before / after table
def before_after():
    rows = []
    for exp in ("E2", "E5"):
        z = L(f"{exp}.npz"); meths = list(z["methods"]); xs = z["x"]
        io, ifx, ior = meths.index("Algo3 M2-original"), meths.index("Algo3 fixed"), meths.index("Ideal gain (oracle)")
        ie6, inn = meths.index("Eq6 (baseline)"), meths.index("No gain correction")
        for tag in ("m2", "paper"):
            chi, rms, lam = z[f"{tag}__chi"], z[f"{tag}__rms"], z[f"{tag}__lam"]
            for j, xv in enumerate(xs):
                snr = lambda i: float(db(np.nanmedian(chi[i, j] / (1 - chi[i, j]))))
                rows.append(dict(experiment=exp, simulator=tag, x=("gamma=%g" % xv) if exp == "E2" else ("K=%d" % xv),
                                 expected_lambda_dB=0.5,
                                 lambda_original_median_dB=round(float(np.nanmedian(lam[0, j])), 3),
                                 lambda_fixed_median_dB=round(float(np.nanmedian(lam[1, j])), 3),
                                 lambda_fixed_bonf_median_dB=round(float(np.nanmedian(lam[2, j])), 3),
                                 abs_err_original_dB=round(abs(float(np.nanmedian(lam[0, j])) - 0.5), 3),
                                 abs_err_fixed_dB=round(abs(float(np.nanmedian(lam[1, j])) - 0.5), 3),
                                 within_0p05_original=round(float(np.mean(np.abs(lam[0, j] - .5) < .05)), 3),
                                 within_0p05_fixed=round(float(np.mean(np.abs(lam[1, j] - .5) < .05)), 3),
                                 within_0p05_fixed_bonf=round(float(np.mean(np.abs(lam[2, j] - .5) < .05)), 3),
                                 rms_gain_err_original_dB=round(float(np.median(rms[io, j])), 4),
                                 rms_gain_err_fixed_dB=round(float(np.median(rms[ifx, j])), 4),
                                 rms_gain_err_none_dB=round(float(np.median(rms[inn, j])), 4),
                                 snr_original_dB=round(snr(io), 2), snr_fixed_dB=round(snr(ifx), 2),
                                 snr_eq6_dB=round(snr(ie6), 2), snr_none_dB=round(snr(inn), 2), snr_ideal_dB=round(snr(ior), 2)))
    with open(os.path.join(R, "before_after.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return rows


if __name__ == "__main__":
    print("lambda example (orig, fixed):", lambda_analysis())
    for r in before_after(): print(r)
    for r in figure6(): print(r)
