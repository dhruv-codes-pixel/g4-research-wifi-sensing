# Git workflow for the M3 gain work (branch `m3-gain-fix`)

Run these from the repository root after copying the delivered files into place.
Each commit holds a real, separable piece of the work, in the order it was done. Review each diff before committing.

```bash
git checkout main && git pull
git checkout -b m3-gain-fix

# 1. Equation -> code mapping (investigation starts here)
git add docs/gain_equation_mapping.md
git commit -m "gain: map Ratnam et al. Eqs. 6-19 and Algorithms 1-3 to the M2 code"

# 2. Simulator extension (csi_sim.py untouched; bit-identical in m2 mode)
git add Code/sim_ext.py
git commit -m "gain: add sim_ext with paper-faithful g1 (<0.1 Hz) and type-(ii) dynamics"

# 3. Instrumented original + corrected Algorithm 3
git add Code/gain_fix.py
git commit -m "gain: fix Algorithm 3 lambda search (0.01 dB grid, Rayleigh acceptance, tie-break, Eq.16 rounding)"

# 4. Unit tests + their results
git add tests/test_gain.py Results/gain/sanity_test.csv Results/gain/unit_test_log.txt
git commit -m "gain: add AGC-only, noise (0.05-0.40 dB) and lambda unit tests"

# 5. 500-trial E2/E5 reruns
git add Code/run_gain_m3.py Results/gain/E2.* Results/gain/E2_* Results/gain/E5.* Results/gain/E5_*
git commit -m "gain: rerun E2 and E5 with 500 trials, seed 42, trial-level npz"

# 6. Figure 6 reproduction
git add Code/digitize_fig6.py Code/make_gain_figures.py Results/gain/FIG6.* Results/gain/FIG6_* \
        Results/gain/fig6_* Results/gain/figure6_comparison.png
git commit -m "gain: digitise paper Fig. 6 and compare with our 500-trial reproduction"

# 7. Before/after + lambda analysis + conclusion
git add Results/gain/before_after.csv Results/gain/lambda_analysis.png docs/gain_report_section.md docs/GIT_COMMITS.md
git commit -m "gain: before/after comparison, lambda analysis and M3 conclusion"

git push -u origin m3-gain-fix
# then open a Pull Request m3-gain-fix -> main on GitHub and request review
```

Note: the `git_commit` field in `E2.json`, `E5.json`, `FIG6.json` reads "unavailable" because the runs were made
outside a git checkout. To record a real hash, rerun inside the branch after commit 3:
`cd Code && python run_gain_m3.py --exp all --trials 500 --seed 42` (~35 min, single core). Results are
bit-identical for the same seed; only the metadata changes. The `code_sha256` field already pins the exact code.
