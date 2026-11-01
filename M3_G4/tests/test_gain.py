"""
tests/test_gain.py -- M3 unit tests for gain correction (Algorithm 3, Ratnam et al. 2024).

Run with pytest (from the repository root):   python -m pytest tests/test_gain.py -q
Run as a script to also write the result tables:  python tests/test_gain.py
    -> Results/gain/sanity_test.csv   (Parts 6-8: AGC-only, noise 0.05-0.40 dB, lambda test)

"noise" = std (dB) of the true CSI power fluctuation Gamma_p added to the observed power series. By the
corrected Lemma 1 (sigma^2 = 18.86 (1-gamma^2)/K): 0.05 dB ~ (K=256, gamma=0.986), 0.10 dB ~ (K=256, gamma=0.94),
0.20 dB ~ (K=64, gamma=0.89), 0.40 dB ~ (K=64, gamma<0.5).  Lambda is the AGC step (true value 0.5 dB).
"""
import os, sys, csv, hashlib
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "Code"))
import csi_sim as c        # noqa: E402
import sim_ext as se       # noqa: E402
import gain_fix as gf      # noqa: E402

M2_CSI_SIM_SHA256 = "9b97503c1e5cc3855b0d0f2b66b5a9c40332f0b68e3010cace776b2ec75d23c6"
LAM, P, N_TRIALS, TOL = 0.5, 300, 100, 0.05
BONF = gf.ALL_FIXES | {"bonf"}
VARIANTS = {"original": dict(variant="original"), "fixed": dict(variant="fixed"), "fixed+bonf": dict(fixes=BONF)}


def synth_power(sigma, rng, drift=False, lam=LAM):
    """Observed power (dB): AGC g2 in {-lam,0,lam} w.p. {.2,.6,.2} + optional slow g1 (<0.1 Hz, std 0.2 dB)
    + Gaussian true-power fluctuation of std sigma."""
    g2 = rng.choice([-lam, 0.0, lam], P, p=[0.2, 0.6, 0.2])
    g1 = np.zeros(P)
    if drift:
        x = se._bandlimited(P, rng, 0.0, 0.1); g1 = 0.2 * x / x.std()
    return g1 + g2 + sigma * rng.standard_normal(P), g1, g2


def run_condition(sigma, drift=False, n=N_TRIALS, seed=2026):
    out = {k: dict(lam=[], rms=[]) for k in VARIANTS}; out["no correction"] = dict(lam=[], rms=[])
    for t in range(n):
        rng = np.random.default_rng([seed, int(sigma * 1000), int(drift), t])
        G, g1, g2 = synth_power(sigma, rng, drift)
        true = g1 + g2
        for k, kw in VARIANTS.items():
            gain, info = gf.gain_algo3_from_power(G, **kw)
            est = np.zeros(P) if gain is None else 20 * np.log10(gain)
            e = est - true; e -= e.mean()
            out[k]["lam"].append(info["lam_hat"]); out[k]["rms"].append(np.sqrt(np.mean(e ** 2)))
        e = -true + true.mean(); out["no correction"]["rms"].append(np.sqrt(np.mean(e ** 2)))
    return {k: dict(lam=np.array(v["lam"]), rms=np.array(v["rms"])) for k, v in out.items()}


# ----------------------------------------------------------------------- integrity of the M2 baseline
def test_csi_sim_unchanged():
    h = hashlib.sha256(open(os.path.join(ROOT, "Code", "csi_sim.py"), "rb").read()).hexdigest()
    assert h == M2_CSI_SIM_SHA256


def test_sim_ext_reproduces_csi_sim():
    for K, g in [(64, 0.9), (256, 0.5), (512, 0.99)]:
        a = c.simulate(K=K, gamma=g, rng=np.random.default_rng(5))
        b = se.simulate(K=K, gamma=g, rng=np.random.default_rng(5))
        assert all(np.array_equal(a[k], b[k]) for k in a if k != "gamma")


def test_building_blocks_match_m2():
    rng = np.random.default_rng(0); x = rng.standard_normal(300)
    assert np.allclose(gf.lpf_ma_fast(x, 60), c.lpf_ma(x, 60), atol=1e-12)
    for v in [0.3, 1.0, 2.5, 7.0, 40.0]:
        assert abs(gf.D(v)[0] - c._D(v)) < 1e-12
    s = c.simulate(K=128, gamma=0.9, rng=rng)
    assert abs(gf.chi_metric_fast(s["ht"], s) - c.chi_metric(s["ht"], s)) < 1e-12


def test_original_variant_is_the_m2_code():
    """gain_fix 'original' == csi_sim.gain_algo3, so every diagnostic of it applies to the M2 code."""
    for i, (K, g) in enumerate([(64, 0.9), (64, 0.99), (256, 0.9), (512, 0.95)] * 2):
        s = c.simulate(K=K, gamma=g, rng=np.random.default_rng(100 + i))
        a = c.gain_algo3(s["ht"]); b, _ = gf.gain_algo3(s["ht"], "original")
        assert np.max(np.abs(a - b)) < 1e-12


# ----------------------------------------------------------------------- diagnosis of the 3.5 dB result
def test_m2_bug_reproduced_and_cause():
    """M2 setting (csi_sim, K=64, gamma=0.9): original picks lambda ~ lambda_max (~3.5 dB); the grid point
    nearest 0.5 dB is rejected by the lambda^2/24 test in every trial."""
    lam_hat, rejected = [], []
    for t in range(20):
        s = c.simulate(K=64, gamma=0.9, rng=np.random.default_rng([42, t]))
        _, info = gf.gain_algo3(s["ht"], "original")
        lam_hat.append(info["lam_hat"]); T = info["table"]
        rejected.append(T[np.argmin(np.abs(T[:, 0] - LAM)), 4] == 0)
    assert 2.5 < np.median(lam_hat) < 4.5 and all(rejected)


def test_original_grid_never_contains_true_step():
    misses = [np.min(np.abs(np.arange(0.05, 1.0001, 0.05) * lm - LAM)) for lm in np.linspace(2.0, 5.0, 301)]
    assert np.median(misses) > 0.02            # typical grid miss of the true step
    G = synth_power(0.0, np.random.default_rng(1))[0]
    assert np.any(np.isclose(gf.gain_algo3_from_power(G)[1]["grid"], LAM))   # FIX-1 grid contains 0.50


def test_eq16_printed_form_is_a_typo():
    """Printed Eq. (16) has '-1/2' outside the floor -> estimates are half-integer multiples of lambda.
    Taken literally, the true lambda always fails the lambda^2/24 test; nearest-integer works."""
    G, _, _ = synth_power(0.0, np.random.default_rng(3))
    _, lit = gf.gain_algo3_from_power(G, variant="original", eq16_literal=True)
    T = lit["table"]; row = T[np.argmin(np.abs(T[:, 0] - LAM))]
    assert row[1] > row[2]                     # mean(Gh^2) > lambda^2/24 -> rejected
    _, ok = gf.gain_algo3_from_power(G, variant="fixed")
    assert abs(ok["lam_hat"] - LAM) < 1e-9


def test_lemma1_constant():
    """Lemma 1 states var(Gamma) = 100(1-g^2)/K; a Taylor expansion of 10log10 gives 18.86(1-g^2)/K."""
    rng = np.random.default_rng(9)
    for K, g in [(64, 0.9), (256, 0.9)]:
        b = np.sqrt(g) * np.exp(1j * rng.uniform(0, 2 * np.pi, K))
        d = np.sqrt((1 - g) / 2) * (rng.standard_normal((20000, K)) + 1j * rng.standard_normal((20000, K)))
        v = np.var(10 * np.log10(np.mean(np.abs(b + d) ** 2, 1)))
        assert abs(v / gf.sigma_gamma_db(K, g) ** 2 - 1) < 0.1
        assert v / gf.sigma_gamma_db(K, g, corrected=False) ** 2 < 0.25


# ----------------------------------------------------------------------- Part 6: AGC only
def test_agc_only_power_series():
    G, _, g2 = synth_power(0.0, np.random.default_rng(11))
    gain, info = gf.gain_algo3_from_power(G, variant="fixed")
    assert abs(info["lam_hat"] - LAM) < TOL
    est = 20 * np.log10(gain); e = est - g2
    assert np.ptp(e) < 1e-9                    # every AGC level recovered exactly (constant offset only)


def test_agc_only_through_csi():
    """Same test through the CSI path (Eq. 7 dB conversion): static channel only, gamma = 1, no d, no g1."""
    rng = np.random.default_rng(12)
    b = c.make_static(64, 1.0, rng); g2 = rng.choice([-0.5, 0, 0.5], P, p=[0.2, 0.6, 0.2])
    ht = (10 ** (g2 / 20))[:, None] * b[None, :]
    gain, info = gf.gain_algo3(ht, "fixed")
    assert abs(info["lam_hat"] - LAM) < TOL and np.ptp(20 * np.log10(gain) - g2) < 1e-9


# ----------------------------------------------------------------------- Parts 7-8: noise and lambda
def test_noise_005_and_010_lambda_recovered():
    for sigma in (0.05, 0.10):
        r = run_condition(sigma, n=60)
        assert abs(np.nanmedian(r["fixed"]["lam"]) - LAM) < TOL
        assert abs(np.nanmedian(r["fixed+bonf"]["lam"]) - LAM) < TOL


def test_noise_020_040_not_worse_than_no_correction():
    """At sigma >= 0.2 dB = 0.4*lambda the AGC level is ambiguous frame by frame (P(error) = 2Q(lambda/2sigma)
    = 21 % at 0.2 dB, 53 % at 0.4 dB): no lambda accuracy is required.  Requirement: a clear gain at 0.2 dB
    and no material harm (<= 0.02 dB RMS) at 0.4 dB.  Observed at 0.4 dB: fixed is ~0.01 dB WORSE than no
    correction (spurious small lambda in a minority of trials); fixed+bonf equals no correction."""
    r = run_condition(0.20, n=60)
    assert np.median(r["fixed"]["rms"]) < np.median(r["no correction"]["rms"]) - 0.05
    r = run_condition(0.40, n=60)
    assert np.median(r["fixed"]["rms"]) <= np.median(r["no correction"]["rms"]) + 0.02
    assert np.median(r["fixed+bonf"]["rms"]) <= np.median(r["no correction"]["rms"]) + 0.005


# ----------------------------------------------------------------------- script mode: write the tables
def write_tables(path):
    rows = []
    for drift in (False, True):
        for sigma in (0.0, 0.05, 0.10, 0.20, 0.40):
            r = run_condition(sigma, drift=drift)
            for k in ("original", "fixed", "fixed+bonf"):
                lam = r[k]["lam"]
                rows.append(dict(test="AGC-only" if sigma == 0 else "noise", large_scale_drift=drift, noise_dB=sigma,
                                 variant=k, true_lambda_dB=LAM, est_lambda_median_dB=round(float(np.nanmedian(lam)), 3),
                                 est_lambda_q25=round(float(np.nanpercentile(lam, 25)), 3),
                                 est_lambda_q75=round(float(np.nanpercentile(lam, 75)), 3),
                                 abs_error_dB=round(abs(float(np.nanmedian(lam)) - LAM), 3),
                                 within_0p05dB_fraction=round(float(np.mean(np.abs(lam - LAM) < TOL)), 3),
                                 rms_gain_error_dB_median=round(float(np.median(r[k]["rms"])), 4),
                                 rms_gain_error_no_correction_dB=round(float(np.median(r["no correction"]["rms"])), 4),
                                 frame_ambiguity_2Q=round(float(2 * __import__("scipy.stats").stats.norm.sf(LAM / (2 * sigma))), 4) if sigma > 0 else 0.0,
                                 trials=N_TRIALS))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return rows


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    ok = 0
    for t in tests:
        try:
            t(); print(f"PASS  {t.__name__}"); ok += 1
        except AssertionError as e:
            print(f"FAIL  {t.__name__}  {e}")
    print(f"{ok}/{len(tests)} tests passed")
    out = os.path.join(ROOT, "Results", "gain", "sanity_test.csv"); os.makedirs(os.path.dirname(out), exist_ok=True)
    rows = write_tables(out)
    print(f"wrote {out}")
    for r in rows:
        print(f"drift={r['large_scale_drift']!s:5s} noise={r['noise_dB']:.2f} {r['variant']:11s} lam_hat={r['est_lambda_median_dB']:.3f} "
              f"[{r['est_lambda_q25']:.2f},{r['est_lambda_q75']:.2f}] err={r['abs_error_dB']:.3f} hit={r['within_0p05dB_fraction']:.2f} "
              f"rms={r['rms_gain_error_dB_median']:.3f} (none {r['rms_gain_error_no_correction_dB']:.3f})")
