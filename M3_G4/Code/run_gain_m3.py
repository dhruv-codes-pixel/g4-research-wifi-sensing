"""
run_gain_m3.py -- M3 gain-correction experiments (E2, E5, Figure 6 reproduction).

Usage (from Code/):
    python run_gain_m3.py --exp E2   --trials 500 --seed 42
    python run_gain_m3.py --exp E5   --trials 500 --seed 42
    python run_gain_m3.py --exp FIG6 --trials 500 --seed 42
    python run_gain_m3.py --exp all  --trials 500 --seed 42

Every trial uses its own generator np.random.default_rng([seed, exp_code, point, trial]); all methods (and
both simulator modes) see the same realisation -> paired, trial-level data for Member E's statistics.
Saved per experiment in ../Results/gain/:  <EXP>.npz (every trial), <EXP>.json (metadata + descriptive
medians/IQR only), <EXP>_log.txt, <EXP>_plot.png.
"""
import argparse, json, hashlib, os, platform, subprocess, sys, time, datetime
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csi_sim as c
import sim_ext as se
import gain_fix as gf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "Results", "gain")
BONF = gf.ALL_FIXES | {"bonf"}
METHODS = ["Eq6 (baseline)", "Algo1", "Algo3 M2-original", "Algo3 fixed", "Algo3 fixed+Bonferroni",
           "No gain correction", "Ideal gain (oracle)"]
EXP_CODE = {"E2": 2, "E5": 5, "FIG6": 6}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unavailable (not run inside a git checkout)"


def one_trial(sim, with_eq8=False):
    f = sim["f"]; ht = sim["ht"]
    hp = ht * np.exp(2j * np.pi * f[None, :] * sim["tau"][:, None]) * np.exp(1j * sim["psi"])[:, None]  # ideal phase
    g_o, i_o = gf.gain_algo3(ht, "original")
    g_f, i_f = gf.gain_algo3(ht, "fixed")
    g_b, i_b = gf.gain_algo3(ht, fixes=BONF)
    gains = {"Eq6 (baseline)": c.gain_eq6(ht), "Algo1": c.gain_algo1(ht), "Algo3 M2-original": g_o,
             "Algo3 fixed": g_f, "Algo3 fixed+Bonferroni": g_b, "No gain correction": np.ones(ht.shape[0]),
             "Ideal gain (oracle)": sim["g"]}
    if with_eq8: gains["Eq8 (baseline)"] = gf.gain_eq8(ht)
    chi = {m: gf.chi_metric_fast(hp / g[:, None], sim) for m, g in gains.items()}
    rms = {m: gf.gain_rms_error_db(g, sim["g"]) for m, g in gains.items()}
    lam = dict(orig=i_o["lam_hat"], fixed=i_f["lam_hat"], bonf=i_b["lam_hat"], lam_max=i_o["lam_max"],
               fb_orig=i_o["fallback"], fb_fixed=i_f["fallback"], fb_bonf=i_b["fallback"])
    return chi, rms, lam


def run(exp, trials, seed, log):
    if exp == "E2":
        xname, xs, fixed, modes, dyns = "gamma", [0.5, 0.7, 0.8, 0.9, 0.95, 0.99], dict(K=64), ["m2", "paper"], ["i"]
    elif exp == "E5":
        xname, xs, fixed, modes, dyns = "K", [64, 128, 256, 512], dict(gamma=0.9), ["m2", "paper"], ["i"]
    else:  # FIG6: paper setting K = 256, paper-faithful g1, gamma grid read off Fig. 6, dynamic types (i), (ii)
        xname, xs, fixed, modes, dyns = "gamma", [0.5, 0.6, 0.7, 0.8, 0.85, 0.87, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99], \
            dict(K=256), ["paper"], ["i", "ii"]
    methods = METHODS + (["Eq8 (baseline)"] if exp == "FIG6" else [])
    data = {}
    t0 = time.time()
    for mode in modes:
        for dyn in dyns:
            tag = f"{mode}_{dyn}" if exp == "FIG6" else mode
            CHI = np.full((len(methods), len(xs), trials), np.nan); RMS = CHI.copy()
            LAM = np.full((3, len(xs), trials), np.nan); LMX = np.full((len(xs), trials), np.nan)
            FB = np.zeros((3, len(xs), trials), bool)
            for j, x in enumerate(xs):
                kw = dict(fixed); kw[xname] = x
                for t in range(trials):
                    rng = np.random.default_rng([seed, EXP_CODE[exp], j, t])
                    sim = se.simulate(rng=rng, g1_mode=mode, dyn=dyn, **kw)
                    chi, rms, lam = one_trial(sim, with_eq8=(exp == "FIG6"))
                    for mi, m in enumerate(methods):
                        CHI[mi, j, t] = chi[m]; RMS[mi, j, t] = rms[m]
                    LAM[:, j, t] = [lam["orig"], lam["fixed"], lam["bonf"]]; LMX[j, t] = lam["lam_max"]
                    FB[:, j, t] = [lam["fb_orig"], lam["fb_fixed"], lam["fb_bonf"]]
                snr = lambda v: db(np.median(v / (1 - v)))
                log(f"[{exp} {tag}] {xname}={x}: " + ", ".join(f"{m}={snr(CHI[mi, j]):.1f}dB" for mi, m in enumerate(methods))
                    + f" | lam_hat median orig={np.nanmedian(LAM[0, j]):.2f} fixed={np.nanmedian(LAM[1, j]):.2f} "
                    f"bonf={np.nanmedian(LAM[2, j]):.2f} | {time.time() - t0:.0f}s")
            data[tag] = dict(chi=CHI, rms=RMS, lam=LAM, lam_max=LMX, fallback=FB)
    return xname, xs, fixed, methods, data


def db(x):
    """10log10 of a (median/quantile of a) linear SNR; non-positive values (chi >= 1, see docs) -> nan."""
    x = np.asarray(x, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(x > 0, 10 * np.log10(np.where(x > 0, x, 1)), np.nan)


def describe(v):
    q = np.nanpercentile(v, [25, 50, 75], axis=-1)
    return dict(q25=q[0].tolist(), median=q[1].tolist(), q75=q[2].tolist())


def save(exp, trials, seed, xname, xs, fixed, methods, data, cmd):
    os.makedirs(OUT, exist_ok=True)
    arrays = {"x": np.array(xs, float), "methods": np.array(methods), "lam_rows": np.array(["orig", "fixed", "fixed+bonf"])}
    for tag, d in data.items():
        for k, v in d.items(): arrays[f"{tag}__{k}"] = v
    np.savez_compressed(os.path.join(OUT, f"{exp}.npz"), **arrays)
    summ = {}
    for tag, d in data.items():
        chi = d["chi"]; s = {}
        for mi, m in enumerate(methods):
            v = chi[mi]
            lin = describe(v / (1 - v))
            s[m] = dict(snr_m2_dB={k: db(q).tolist() for k, q in lin.items()},
                        snr_paper_linear=describe(v ** 2 / (1 - v ** 2)),
                        rms_gain_err_dB=describe(d["rms"][mi]))
        s["lambda_hat"] = {r: describe(d["lam"][i]) for i, r in enumerate(["orig", "fixed", "fixed+bonf"])}
        s["lambda_within_0.05dB_of_0.5_fraction"] = {r: np.mean(np.abs(d["lam"][i] - 0.5) < 0.05, axis=-1).tolist()
                                                    for i, r in enumerate(["orig", "fixed", "fixed+bonf"])}
        s["fallback_fraction"] = {r: d["fallback"][i].mean(-1).tolist() for i, r in enumerate(["orig", "fixed", "fixed+bonf"])}
        summ[tag] = s
    meta = dict(experiment=exp, command=cmd, seed=seed, trials=trials, x_name=xname, x_values=xs, fixed_params=fixed,
                sim_defaults=dict(P=300, Trep_s=c.TREP, Ts_s=c.TS, agc_dB=[-0.5, 0, 0.5], agc_p=[0.2, 0.6, 0.2],
                                  g1_std_dB=0.2, tau_max_s=1e-7, static="csi_sim.make_static (6-tap, 200 ns)"),
                simulator_modes={"m2": "csi_sim.simulate (unchanged M2 simulator; g1 = 21-tap MA, ~0.5 Hz)",
                                 "paper": "sim_ext g1 band-limited to [0, 0.1] Hz (paper Sec. V-B), same std normalisation"},
                dynamic_types={"i": "d ~ CN(0,1-gamma) iid", "ii": "single moving path, Doppler support [0.5,1] Hz"},
                per_trial_rng="np.random.default_rng([seed, exp_code, point_index, trial])",
                exp_code=EXP_CODE[exp], methods=methods,
                metric=dict(chi="Eq.(29) as in csi_sim.chi_metric", snr_m2="chi/(1-chi) (M2 convention)",
                            snr_paper="chi^2/(1-chi^2) (paper Fig. 6/7 axis)", rms_gain_err_dB="RMS of 20log10(g_hat/g), offset removed"),
                ideal_phase_correction=True,
                npz_layout="<tag>__chi / __rms : [method, point, trial];  <tag>__lam / __fallback : [orig|fixed|fixed+bonf, point, trial]",
                git_commit=git_commit(),
                code_sha256={f: sha(os.path.join(HERE, f)) for f in ["csi_sim.py", "sim_ext.py", "gain_fix.py", "run_gain_m3.py"]},
                versions=dict(python=platform.python_version(), numpy=np.__version__),
                created_utc=datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
                summary_descriptive_only=summ,
                note="Descriptive medians/IQR only. Bootstrap CIs and Wilcoxon tests are computed by Member E from the .npz.")
    json.dump(meta, open(os.path.join(OUT, f"{exp}.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


COL = {"Eq6 (baseline)": "tab:blue", "Algo1": "tab:orange", "Algo3 M2-original": "tab:green", "Algo3 fixed": "tab:purple",
       "Algo3 fixed+Bonferroni": "tab:brown", "No gain correction": "tab:red", "Ideal gain (oracle)": "black",
       "Eq8 (baseline)": "tab:olive"}


def plot_e(exp, xname, xs, methods, data):
    tags = list(data); fig, ax = plt.subplots(2, len(tags), figsize=(6.2 * len(tags), 8), squeeze=False)
    names = {"m2": "M2 simulator (csi_sim, g1 ~0.5 Hz)", "paper": "paper-faithful g1 (<0.1 Hz)"}
    for ci, tag in enumerate(tags):
        d = data[tag]
        for mi, m in enumerate(methods):
            if m == "Algo3 fixed+Bonferroni": continue
            v = d["chi"][mi] / (1 - d["chi"][mi])
            q1, md, q3 = db(np.nanpercentile(v, [25, 50, 75], axis=1))
            ls = "--" if m == "Ideal gain (oracle)" else "-"
            ax[0, ci].plot(xs, md, ls, marker="o", ms=4, color=COL[m], label=m); ax[0, ci].fill_between(xs, q1, q3, color=COL[m], alpha=0.12)
            r = d["rms"][mi]
            if m != "Ideal gain (oracle)":
                ax[1, ci].plot(xs, np.nanmedian(r, 1), ls, marker="o", ms=4, color=COL[m], label=m)
        ax[0, ci].set_title(f"{exp}: {names.get(tag, tag)}", fontsize=10)
        ax[0, ci].set_ylabel("Median SNR chi/(1-chi) (dB)"); ax[1, ci].set_ylabel("Median RMS gain error (dB)")
        for a in ax[:, ci]:
            a.grid(alpha=.3); a.set_xlabel(xname if xname == "gamma" else "subcarriers K")
            if xname == "K": a.set_xscale("log", base=2); a.set_xticks(xs); a.set_xticklabels([str(x) for x in xs])
        ax[0, ci].legend(fontsize=7); ax[1, ci].legend(fontsize=7)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, f"{exp}_plot.png"), dpi=150); plt.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", default="all", choices=["E2", "E5", "FIG6", "all"])
    ap.add_argument("--trials", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    cmd = " ".join(["python"] + sys.argv)
    for exp in (["E2", "E5", "FIG6"] if a.exp == "all" else [a.exp]):
        os.makedirs(OUT, exist_ok=True)
        lf = open(os.path.join(OUT, f"{exp}_log.txt"), "w")
        def log(s, lf=lf):
            line = f"{datetime.datetime.utcnow().isoformat(timespec='seconds')}Z  {s}"; print(line, flush=True); lf.write(line + "\n"); lf.flush()
        log(f"command: {cmd}"); log(f"seed={a.seed} trials={a.trials} csi_sim.py sha256={sha(os.path.join(HERE, 'csi_sim.py'))}")
        xname, xs, fixed, methods, data = run(exp, a.trials, a.seed, log)
        save(exp, a.trials, a.seed, xname, xs, fixed, methods, data, cmd)
        if exp != "FIG6": plot_e(exp, xname, xs, methods, data)
        log("done; saved npz/json" + ("" if exp == "FIG6" else "/png")); lf.close()
