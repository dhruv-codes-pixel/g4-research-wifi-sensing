# M3 – Gain correction (Algorithm 3): equation → code mapping and investigation

Base paper: V. V. Ratnam et al., "Optimal preprocessing of WiFi CSI for sensing applications", IEEE TWC 23(9), 2024
(arXiv:2307.12126v2 used for equations and figures). Scope: Sec. III, Eqs. 6–19, Algorithms 1–3, Fig. 6, Sec. V-B.

Files: M2 code `Code/csi_sim.py` (unchanged, sha256 `9b97503c…d23c6`); M3 code `Code/gain_fix.py`, `Code/sim_ext.py`,
`Code/run_gain_m3.py`, `Code/digitize_fig6.py`, `Code/make_gain_figures.py`; tests `tests/test_gain.py`.

## 1. Equation → code mapping

✓ = M2 code matches the paper; ✗ = mismatch found; — = not in M2.

| Paper | Meaning | M2 code (`csi_sim.py`) | M3 code | Status |
|---|---|---|---|---|
| Eq. 2 | Observed CSI h̃ = g·h·e^{−j2πf_kτ}·e^{−jψ} | `simulate` L43 | `sim_ext.simulate` | ✓ |
| Eq. 3 / Eq. 4 | Corrected CSI, h̄ = h̃ / ĝ | `run_gain_trial` L171–173 (`hp / gain`) | `run_gain_m3.one_trial` | ✓ |
| Eq. 5 | 20·log10 g = g⁽¹⁾ + g⁽²⁾ (dB) | `simulate` L40 | same | ✓ |
| Sec. V-B, g⁽¹⁾ | Gaussian process, std 0.2 dB, Doppler support [0, 0.1] Hz | `simulate` L37–38: 21-tap moving average at T_rep = 0.1 s | `sim_ext` `g1_mode="paper"` (FFT band-limit to 0.1 Hz, same std normalisation) | ✗ about 60 % of g⁽¹⁾ power lies above 0.1 Hz (first spectral null ≈ 0.48 Hz) |
| Sec. V-B, g⁽²⁾ | AGC ∈ {−0.5, 0, 0.5} dB w.p. {0.2, 0.6, 0.2} | `simulate` L39 | same | ✓ |
| Eq. 6 | Baseline: ĝ = √(mean_k abs(h̃)²) | `gain_eq6` L47–48 | reused | ✓ |
| Eq. 7 | Γ̃_p = 10·log10 mean_k abs(h̃)² | L59 (Algo 1), L78 (Algo 3) | `gain_fix.csi_power_db` | ✓ |
| Eq. 8 | Baseline: DBSCAN on Γ̃ (ε = 0.15, min-points 1), cluster mean | — | `gain_fix.gain_eq8` | — added (needed for Fig. 6) |
| Eq. 9 | True power Γ_p | (tests only) | `tests/test_gain.py` | ✓ |
| Lemma 1 | var Γ_p = 100(1−γ²)/K | not used by code | `gain_fix.sigma_gamma_db` | ✗ in the paper: correct constant is (10/ln10)² = 18.86 (test `test_lemma1_constant`) |
| Eq. 10 | ΔΓ̃_p = Γ̃_p − Γ̃_{p−1} | L60 | reused | ✓ |
| Eq. 12 | ĝ⁽²⁾ by cumulative sum of DBSCAN cluster means | L61–65 | reused | ✓ |
| Eq. 13 | ĝ⁽¹⁾ = LPF_{0.1 Hz}{Γ̃ − ĝ⁽²⁾} | L66 (`lpf_ma`, one-sided 60) | reused | ✓ literal; see §3-E on the cutoff |
| Algorithm 1 | Gain estimation, arbitrary AGC grid | `gain_algo1` L57–67 (ε = 0.2, as in the Algorithm 1 box) | reused unchanged | ✓ |
| Lemma 2 / Eq. 14 | Ξ̃_p = exp(j2πΓ̃_p/λ) | L82 | `gain_fix.algo2_g1` | ✓ |
| Eq. 15 | Ξ̃ rotates slowly through g⁽¹⁾ only | (rationale) | — | ✓ |
| **Algorithm 2** | ĝ⁽¹⁾ = U[∠ LPF{Ξ̃}]·λ/2π; moving average, one-sided width 6/T_rep | inline L83–84 | `gain_fix.algo2_g1` (vectorised, identical) | ✓ (see §2) |
| Eq. 16 | ML AGC estimate. Printed: λ(⌊x/λ + ½⌋ − ½) | L85 `lam*np.round(x/lam)` | `floor(x/λ + ½)` (FIX-4) | ✓ M2 is right; the printed "−½" is a typo (`test_eq16_printed_form_is_a_typo`) |
| Eq. 17 | Γ̂ = Γ̃ − ĝ⁽¹⁾ − ĝ⁽²⁾ | L86 | same | ✓ |
| Lemma 3 / Eq. 18 | σ̂²_Γ = −λ²/(2π²)·ln abs(mean exp(j2πΓ̂/λ)) | L89–90 | same | ✓ |
| Lemma 4 / Eq. 19 | Distortion λ²·D(λ/σ_Γ), D(x) = Σ_z [Q((z−½)x) − Q((z+½)x)]·z² | `_D` L71–73 (abs z ≤ 30) | `gain_fix.D` (vectorised, identical) | ✓ |
| Algorithm 3, line 2 | λ_max = 1.5(max Γ̃ − min Γ̃) | L79 | same | ✓ |
| Algorithm 3, loop | λ ∈ (0.05 : 0.05 : 1)·λ_max | L81 | 0.01 dB multiples in the same range (FIX-1) | ✓ literal, but too coarse (§3-C/D) |
| Algorithm 3, test | accept λ only if Σ_p Γ̂²/P ≤ λ²/24 | L87–88 | Rayleigh test P·m² > ln(1/α), α = 10⁻³ (FIX-2) | ✓ literal; this is the **root cause** (§3-B) |
| Algorithm 3, objective | Obj(λ) = σ̂²_Γ + λ²·D(λ/σ̂_Γ) | L91 | same | ✓ |
| Algorithm 3, selection | λ̂ = argmin Obj | L92–93 (first minimum) | ties → largest λ (FIX-3) | ✓ |
| Algorithm 3, output | ĝ = 10^{(ĝ⁽¹⁾+ĝ⁽²⁾)/20} | L93 | same | ✓ |
| (not in paper) | No λ accepted | L94: silently returns Eq. 6 | same fallback, but flagged in `info["fallback"]` (FIX-5) | ✗ undocumented in M2 |
| Eq. 29 | χ (squared correlation of the dynamic part) | `chi_metric` L141–151 | `gain_fix.chi_metric_fast` (identical, cached) | ✓ |
| Fig. 6 y-axis | SNR = χ²/(1−χ²) | M2 uses χ/(1−χ) (L153–154) | both saved | ✗ convention differs (M2 kept for continuity) |

The `variant="original"` path of `gain_fix.gain_algo3` re-implements L75–94 and reproduces `csi_sim.gain_algo3` to
2×10⁻¹⁶ (`test_original_variant_is_the_m2_code`). Every diagnostic below therefore applies to the M2 code itself.

## 2. Is Algorithm 2 required?

Yes. Algorithm 2 is a step inside every λ iteration of Algorithm 3, because ĝ⁽¹⁾ must be known before Eq. 16 can
round the residual to the AGC grid. M2 already implements it inline (L82–84), and it matches the paper: Ξ̃ from
Eq. 14, a moving average with one-sided width 6/T_rep = 60 frames applied to the real and imaginary parts, then
the unwrapped angle times λ/2π. No new implementation was needed. M3 only factors it into `gain_fix.algo2_g1`
(vectorised, numerically identical) so it can be tested and instrumented. Its accuracy limits Algorithm 3 (§3-E).

## 3. Investigation of λ̂ ≈ 3.5 dB instead of 0.5 dB

The bug reproduces on the unchanged M2 simulator at K=64, γ=0.9: the median λ̂ is 3.54 dB, i.e. λ̂ ≈ λ_max.
In 100 % of trials the grid point nearest 0.5 dB is rejected (`test_m2_bug_reproduced_and_cause`).
`lambda_analysis.png`, top-left, shows one trial.

**A. Objective function.** Eq. 18, Eq. 19 and Obj(λ) match the paper term by term, including signs, normalisation,
natural log and the abs z ≤ 30 truncation (`test_building_blocks_match_m2`). The objective is never evaluated near
0.5 dB because those λ are filtered out first (B). There is one property of the objective worth knowing. The λ² in
Eq. 18 makes σ̂² small for small λ even when m is low, so the objective is biased towards small λ whenever m is near
its noise floor. This is the source of the occasional sub-multiple choices after the fix.

**B. λ estimation: acceptance test (root cause).** Algorithm 3 keeps λ only if mean(Γ̂²) ≤ λ²/24. The value λ²/24
is half the variance of a uniform distribution on [−λ/2, λ/2]. For λ = 0.5 dB, with perfect g⁽¹⁾, the test passes
only if the true power fluctuation σ_Γ < 0.105 dB. Using the corrected Lemma 1, σ_Γ is:

| σ_Γ (dB) | γ = 0.5 | 0.7 | 0.8 | 0.9 | 0.95 | 0.99 |
|---|---|---|---|---|---|---|
| K = 64 | 0.470 | 0.388 | 0.326 | 0.237 | 0.170 | 0.077 ✓ |
| K = 128 | 0.332 | 0.274 | 0.230 | 0.167 | 0.120 ⊕ | 0.054 ✓ |
| K = 256 | 0.235 | 0.194 | 0.163 | 0.118 ⊕ | 0.085 ✓ | 0.038 ✓ |
| K = 512 | 0.166 | 0.137 ⊕ | 0.115 ⊕ | 0.084 ✓ | 0.060 ✓ | 0.027 ✓ |

✓ passes the paper's test (< 0.105 dB). ⊕ passes only the FIX-2 Rayleigh test (< 0.155 dB).

At M2's operating point (K=64, γ=0.9) the true λ cannot pass. Every small λ is rejected; only λ close to λ_max
(≈ 1.5 × the 2–3 dB range of Γ̃ ≈ 3.5 dB) passes. At that λ, Eq. 16 gives ĝ⁽²⁾ = 0, so Algorithm 3 degenerates to
"g⁽¹⁾ only". This is exactly why M2's E2/E5 show Algorithm 3 ≈ No gain correction. The test is also more
conservative than needed. The circular statistic of Eq. 18 estimates σ_Γ correctly even when Γ̂ wraps, as long as
m stays well above its noise level 1/√P ≈ 0.06. FIX-2 therefore tests that statistic directly: a Rayleigh test of
non-uniformity, P·m² > ln(10³), i.e. m > 0.152 for P = 300.

**C. Search range.** The range [0.05, 1]·λ_max is exactly the paper's, and it contains 0.5 dB. It is not the cause.

**D. Search resolution.** The paper's 20-point grid has a step of 0.05·λ_max ≈ 0.14–0.21 dB, so 0.5 dB is never on
it; the typical miss is > 0.02 dB. A wrong λ (e.g. 0.585) shifts every AGC frame by ∓0.085 dB, which inflates
mean(Γ̂²) and makes rejection even more certain. FIX-1 uses multiples of 0.01 dB. Ablation (40 trials, K=256, γ=0.9):
FIX-1 alone changes nothing (λ̂ = 2.94 dB). FIX-2 alone gives 0.57 dB. FIX-1 + FIX-2 gives 0.47 dB.

**E. Low-pass filter (Algorithm 2).** It is literal in M2: a moving average with one-sided width 6/T_rep = 60
frames (121 taps, 12.1 s), applied to Ξ̃, not to Γ̃, before the rounding of Eq. 16. Three findings:

1. That filter's −3 dB cutoff is ≈ 0.037 Hz, not the "0.1 Hz" in its name.
2. M2's g⁽¹⁾ is about 5× faster than the paper's model (table row "Sec. V-B, g⁽¹⁾"). Algorithm 2 cannot track it,
   and the residual inflates Γ̂.
3. A filter with a true 0.1 Hz cutoff (one-sided width 22) was tested. It tracks g⁽¹⁾ better, but it makes
   Algorithm 2 follow the noise at small λ, and λ̂ collapses to about 0.16–0.21 dB in ≥ 87 % of trials.

The paper's width was therefore kept. Even with the paper-faithful g⁽¹⁾, Algorithm 2's tracking error lowers m at
λ = 0.5 from the ideal 0.33 to about 0.22 (K=256, γ=0.9). This is the main remaining limitation.

**F. Rounding.** The order matches the paper: Algorithm 2 (g⁽¹⁾), then round (Eq. 16), then residual (Eq. 17), then
test and objective, then argmin. Rounding happens inside each λ hypothesis, never after optimisation.
`np.round` (half-to-even) and the paper's ⌊x + ½⌋ (half-up) differ only at exact ties; FIX-4 follows the paper.
The printed Eq. 16 has "− ½" outside the floor, which would give half-integer multiples of λ (contradicting
g⁽²⁾ = zλ). Taken literally, the true λ always fails the test, so it is a typo.

**Other checks.**
- Lemma 1's constant is wrong by (ln 10)² ≈ 5.3. The simulated variance matches 18.86(1−γ²)/K to < 10 %.
  The algorithm is unaffected (σ̂_Γ comes from data via Eq. 18), but the paper's theory overstates σ_Γ by 2.3×.
- Tried and rejected: (i) a leave-one-out moving average in Algorithm 2, which brought no improvement in λ̂ or SNR;
  (ii) Bonferroni correction of FIX-2 over the about 280 grid points. Bonferroni removes spurious sub-multiple λ but
  rejects the true λ more often, which lowers downstream SNR. It is kept as the optional variant
  `Algo3 fixed+Bonferroni` and reported alongside.

## 4. Changes in `gain_fix.py` (original → fixed)

| ID | Original (paper / M2) | Fixed | Why |
|---|---|---|---|
| FIX-1 | λ grid (0.05:0.05:1)·λ_max, 20 points | multiples of 0.01 dB in [0.05, 1]·λ_max | 0.5 dB never on the grid |
| FIX-2 | accept if mean(Γ̂²) ≤ λ²/24 | accept if P·m² > ln(1/α), α = 10⁻³ | λ²/24 rejects the true λ whenever σ_Γ > 0.105 dB, although Eq. 18 is still valid |
| FIX-3 | first minimum of Obj | among ties (≤ 10⁻⁹ dB²), the largest λ | λ and λ/n correct the AGC identically without noise |
| FIX-4 | `np.round` | ⌊x + ½⌋ | as printed in Eq. 16 (minus the typo) |
| FIX-5 | silent Eq. 6 fallback | fallback flagged in `info` | the paper does not define this case |

Algorithm 2, Eq. 14, Eq. 17, Eq. 18, Eq. 19 and the objective are unchanged.

## 5. Results

### 5.1 Unit tests (`tests/test_gain.py`, 12/12 pass; table `Results/gain/sanity_test.csv`)

Synthetic observed power = AGC (λ = 0.5 dB) + Gaussian true-power fluctuation ("noise"), P = 300, 100 trials per row.
Values are the median λ̂ in dB, with the fraction of trials within ±0.05 dB of 0.5 in brackets.

| Condition | True λ | Original | Fixed | Fixed + Bonferroni | Abs. error, fixed | Frame ambiguity 2Q(λ/2σ) |
|---|---|---|---|---|---|---|
| AGC only (clean) | 0.5 | 0.525 (1.00) | **0.500** (1.00) | 0.500 (1.00) | 0.000 | 0 |
| noise 0.05 dB | 0.5 | 0.474 (0.95) | **0.500** (0.91) | 0.500 (1.00) | 0.000 | 0.0000 |
| noise 0.10 dB | 0.5 | 0.544 (0.59) | **0.510** (0.77) | 0.520 (0.99) | 0.010 | 0.012 |
| noise 0.20 dB | 0.5 | 2.874 (0.00) | 0.600 (0.05) | 2.840 (0.03) | 0.100 | 0.21 |
| noise 0.40 dB | 0.5 | 4.176 (0.00) | 4.040 (0.01) | 4.170 (0.00) | 3.54 | 0.53 |

- **AGC only:** the fixed version recovers every AGC level exactly, through both the power path and the full CSI
  path (`test_agc_only_*`).
- **0.05–0.10 dB:** λ is recovered within the 0.05 dB tolerance.
- **0.20–0.40 dB:** the AGC level is ambiguous frame by frame, so no λ accuracy is possible.
  - At 0.20 dB the fixed version still cuts the RMS gain error from 0.316 to 0.226 dB.
  - At 0.40 dB it is 0.006 dB worse than no correction; the Bonferroni variant equals no correction.
- **With slow drift added (g⁽¹⁾, < 0.1 Hz):** spurious sub-multiples appear more often (median 0.25 dB at
  noise 0.05; 0.465 dB at noise 0.10). The RMS gain error still falls from about 0.34 to 0.08–0.13 dB, because
  λ/2 corrects the AGC almost as well as λ.

### 5.2 Before / after (`Results/gain/before_after.csv`)

| Setting | Expected λ | Original λ̂ | Fixed λ̂ (within ±0.05) | RMS gain error, orig → fixed | Original SNR | Fixed SNR | Ideal-gain SNR |
|---|---|---|---|---|---|---|---|
| M2 bug case: K=64, γ=0.9, M2 simulator | 0.5 | 3.59 | 3.42 (3 %); σ_Γ = 0.237 dB, not identifiable | 0.366 → 0.354 dB | 16.6 dB | 17.3 dB | 20.0 dB |
| K=64, γ=0.99, M2 simulator | 0.5 | 2.92 | **0.47** (33 %) | 0.365 → 0.132 dB | 7.4 dB | **15.6 dB** | 19.5 dB |
| K=512, γ=0.9, M2 simulator | 0.5 | 2.93 | **0.48** (31 %) | 0.364 → 0.139 dB | 16.7 dB | **21.7 dB** | 23.9 dB |
| K=256, γ=0.9, paper g⁽¹⁾ (paper's setting) | 0.5 | 2.96 | **0.48** (38 %) | 0.342 → 0.149 dB | 17.2 dB | **21.3 dB** | 22.9 dB |
| K=512, γ=0.9, paper g⁽¹⁾ | 0.5 | 2.81 | **0.46** (49 %) | 0.346 → 0.114 dB | 17.1 dB | **22.6 dB** | 23.9 dB |

The median λ̂ moves from about 3 dB to about 0.5 dB, but the per-trial spread is wide: only 31–49 % of trials fall
within ±0.05 dB. The remaining trials are mostly spurious sub-multiples (0.2–0.35 dB), which still remove most of
the AGC error. At γ ≤ 0.8 and K = 64 the fixed version is 0.0–0.2 dB below the original (both are ≈ no correction).

SNR = median χ/(1−χ) (M2 convention), 500 trials, seed 42.

### 5.3 E2 and E5 (500 trials, seed 42, `E2_plot.png`, `E5_plot.png`)

- **E2 (K = 64).** The fixed version changes nothing for γ ≤ 0.8, because σ_Γ ≥ 0.33 dB and the AGC is not
  identifiable. It gains about 0.7 dB at γ = 0.9 and 1.0–2.0 dB at γ = 0.95. At γ = 0.99 it gains +8.2 dB on the
  M2 simulator and +9.2 dB with the paper g⁽¹⁾. The original stays within 0.5 dB of no correction at every γ.
- **E5 (γ = 0.9).** This is the clearest result. The original is flat at about 16.8 dB for all K, which is the
  anomaly flagged in M2. The fixed version now improves with K as theory predicts (σ_Γ ∝ 1/√K):
  - M2 simulator: 17.3 → 17.7 → 19.6 → 21.7 dB;
  - paper g⁽¹⁾: 17.4 → 18.6 → 21.3 → 22.6 dB.
  Its median λ̂ reaches 0.46–0.48 dB at K ≥ 256.
- **Eq. 6 remains best in our metric** (19.2 → 23.4 dB). The reason is that Eq. 6's only error is to remove the
  true power fluctuation. Its RMS gain error equals σ_Γ exactly (0.235/0.167/0.118/0.084 dB against the corrected
  Lemma 1's 0.237/0.167/0.118/0.084 dB), and that error costs little under Eq. 29 with P = 300. Even ideal gain
  correction is only 0.5–1.0 dB above Eq. 6 (finite-P ceiling of Eq. 29; see §5.4).

### 5.4 Figure 6 (`figure6_comparison.png`, `fig6_table.csv`, `fig6_paper_digitised.csv`, `FIG6.npz`)

**Setup.** This matches the paper's Sec. V-B setting: K = 256, P = 300, T_rep = 0.1 s, paper-faithful g⁽¹⁾,
AGC {−0.5, 0, 0.5} dB, τ ~ U[0, 100 ns), ideal phase correction. It uses the paper's 13 γ values, 500 trials per
point, seed 42, and both dynamic types. The paper's curves were digitised from the arXiv raster images with
`digitize_fig6.py` (98 points, ±3 %; overlays saved). Values below are median χ²/(1−χ²), the paper's axis.

| Type, γ | Paper Eq. 6 | Paper Algo 3 | Ours Eq. 6 | Ours Algo 3 fixed | Ours Algo 3 original | Ours ideal gain |
|---|---|---|---|---|---|---|
| (i), 0.6 | 590 | 689 | 85 | 71 | 65 | 100 |
| (i), 0.8 | 518 | 1236 | 87 | 61 | 44 | 107 |
| (i), 0.9 | 494 | 1218 | 84 | 71 | 26 | 102 |
| (i), 0.96 | 478 | 491 | 87 | 69 | 11 | 102 |
| (ii), 0.6 | 6.7 | 35 | 6.3 | 6.0 | 40 | 5857 |
| (ii), 0.8 | 8.2 | 23 | 12 | 12 | 41 | 6387 |
| (ii), 0.9 | 9.3 | 16 | 24 | 24 | 27 | 6193 |

**Type (i): where it matches.**
- Ordering for γ ≤ 0.92: Eq. 6 above Algorithm 1 and Eq. 8, which fall with γ.
- The fixed Algorithm 3 stays close to Eq. 6 and well above Algorithm 1/Eq. 8 at γ ≥ 0.8.
  - γ = 0.9: 71 against 24/23. In the paper, 1218 against 6/5.
- The fixed version is the only Algorithm 3 that stays useful above γ ≈ 0.9, as in the paper.
- The original collapses with γ: 26 at γ = 0.9, 11 at γ = 0.96.

**Type (i): where it differs.**
1. **Absolute scale.** Our numbers are 6–17× lower. The cause is Eq. 29 as printed with P = 300: even ideal gain
   correction peaks at about 100, below the paper's Eq. 6 (≈ 500) and Algorithm 3 (≈ 1300). The static part b̂ is
   a sample mean over P frames, which caps χ at about 1 − 1/P, and the denominator uses the expected rather than the
   sample dynamic energy. The paper's numbers therefore cannot come from Eq. 29 exactly as printed. Its
   implementation must differ in a way the paper does not report.
2. **Algorithm 3 above Eq. 6** (paper: up to 2.5×, +4 dB) is not reproduced. Our fixed Algorithm 3 stays 0.4–1.9 dB
   below Eq. 6, and the ideal gain is only 0.6–1.0 dB above Eq. 6. Eq. 6's only error in type (i) is the true power
   fluctuation σ_Γ (0.118 dB at K = 256, γ = 0.9). This costs little under this metric, so no gain method can beat
   Eq. 6 by much here.
3. **Algorithm 1 and Eq. 8 recovery at γ ≥ 0.94** in the paper (to about 400) is not seen; ours stay at
   no-correction level. With 0.2 dB of g⁽¹⁾ drift, the AGC levels of Γ̃ smear into a continuum and single-linkage
   DBSCAN (min-points 1) chains them into one cluster. Our Eq. 8 equals no correction at every γ.

**Type (ii).**
- Eq. 6 is the worst method at low γ, as in the paper (6.3 against 6.7 at γ = 0.6), because the sensing signal sits
  in the CSI power.
- Our **original** Algorithm 3 has the paper's Algorithm 3 magnitude at low γ (40 against 41 at γ = 0.5). It
  effectively applies no AGC correction, and that preserves the power signal.
- The **fixed** Algorithm 3 drops to the Eq. 6 level for γ ≤ 0.87. Type-(ii) power fluctuations are large
  (several dB) and correlated over about 5–10 frames (Doppler 0.5–1 Hz). That violates the independent-frames
  assumption behind FIX-2's Rayleigh threshold, as well as Lemma 1. Spurious small λ are accepted, and the
  "AGC correction" then removes the sensing signal.
- **FIX-2 is therefore valid for type (i) only.** A next step is an effective-sample-size correction of the
  threshold (ln(1/α)/N_eff, with N_eff from the autocorrelation of Γ̃).
- The ideal gain gives about 6000 because a rank-1 dynamic component is almost perfectly recovered once the gain is
  exact.

## 6. Conclusion: **Partially reproduced**

**Reproduced and fixed.**
- The ~3.5 dB result is explained. Algorithm 3's acceptance test (Σ Γ̂²/P ≤ λ²/24) rejects the true 0.5 dB step
  whenever the true power fluctuation exceeds ~0.105 dB, and the 20-point grid never contains 0.5 dB. Only λ ≈ λ_max
  survives, which applies no AGC correction. That is why M2 found Algorithm 3 ≈ "No gain correction".
- With FIX-1 to FIX-5, in the regime the algorithm is designed for (i.i.d. dynamics, σ_Γ ≲ 0.15 dB: K ≥ 256 at
  γ = 0.9, or K = 64 at γ = 0.99):
  - median λ̂ = 0.46–0.48 dB (original 2.8–3.6 dB);
  - RMS gain error 0.34–0.37 → 0.11–0.18 dB;
  - downstream SNR +4.1 to +9.2 dB over the original, within 1.3–3.9 dB of ideal gain correction.
- Gain-correction quality now grows with bandwidth K, which M2 expected but did not observe.
- All unit tests pass, including exact recovery in the AGC-only case.

**Not reproduced.**
1. Fig. 6 magnitudes and "Algorithm 3 > Eq. 6" in type (i). This is not possible under Eq. 29 as printed, because
   the finite-P ceiling sits below the paper's values.
2. Exact λ in every trial. Per-trial λ̂ is broad (31–49 % within ±0.05 dB). Algorithm 2's g⁽¹⁾ tracking error and
   the λ²-bias of Eq. 18 favour sub-multiples.
3. Type-(ii) behaviour of the fixed version (see §5.4).
4. Algorithm 1 and Eq. 8 recovery at γ ≥ 0.94.

Algorithm 3 is fundamentally limited when σ_Γ ≥ 0.4λ (frame-level ambiguity of 21–53 %). This covers M2's original
E2 operating point (K = 64, γ ≤ 0.9).

## 7. Limitations

- Simulated data only. The static channel is the M2 6-tap model, not 802.11ax Model-C.
- The paper's g⁽¹⁾ spectrum shape inside [0, 0.1] Hz is unspecified; a flat spectrum was assumed (worst case for
  Algorithm 2).
- Eq. 29 is used exactly as in M2. The paper's exact implementation of its metric is unknown.
- The arXiv v2 text was used. Eq. 16, Algorithm 3 and Lemma 1 should be re-checked against the IEEE published PDF.
- FIX-2 assumes independent frames (§5.4). Bonferroni is optional.
- Only descriptive statistics (median/IQR) are reported here. Inference is Member E's task.

## 8. Handover to Member E (statistics)

- **Files:** `Results/gain/E2.npz`, `E5.npz`, `FIG6.npz` (trial level) and `E2.json`, `E5.json`, `FIG6.json`
  (seed, trials, parameters, code sha256, method order, layout).
- **Array layout:**
  - `<tag>__chi` and `<tag>__rms` have shape [method, point, trial].
  - `<tag>__lam` and `<tag>__fallback` have shape [original | fixed | fixed+Bonferroni, point, trial].
  - Tags: `m2`, `paper` (E2/E5); `paper_i`, `paper_ii` (FIG6). `x` holds the γ or K values; `methods` holds the
    method names.
- **Pairing.** All methods (and both simulator tags) share the same realisation per trial
  (`default_rng([42, exp, point, trial])`). Use paired tests: Wilcoxon signed-rank on per-trial differences.
- **Main comparisons:**
  - "Algo3 fixed" vs "Algo3 M2-original";
  - "Algo3 fixed" vs "Eq6 (baseline)";
  - "Algo3 fixed" vs "No gain correction".
  Run them per point, on (a) χ and (b) RMS gain error.
- **Caution:** with ideal or near-ideal correction, χ can exceed 1 (Eq. 29 normalisation), which makes χ/(1−χ)
  negative. Test on χ itself or on RMS gain error, not on the dB SNR. Bootstrap medians of χ/(1−χ) are fine if
  computed on the linear values.

## 9. Material for the final report (Part 19)

**Method.** Algorithm 3 estimates the receiver gain of each CSI frame as a slow large-scale part g⁽¹⁾ plus an AGC
part g⁽²⁾ on an unknown uniform grid of step λ. For each candidate λ it does the following: phase-wraps the
observed power (Eq. 14); tracks g⁽¹⁾ by low-pass filtering (Algorithm 2); rounds the residual to the grid (Eq. 16);
estimates the residual spread (Eq. 18); and scores the candidate by spread plus quantisation distortion (Eq. 19).
The M2 implementation followed the paper literally but returned λ ≈ 3.5 dB instead of 0.5 dB. We traced this to
the λ²/24 acceptance test and the 20-point λ grid. We replaced them with a Rayleigh non-uniformity test on the
Eq. 18 statistic and a 0.01 dB grid, plus a tie-break, the paper's rounding rule and an explicit fallback flag.

**Setup.** P = 300 frames, T_rep = 0.1 s, AGC {−0.5, 0, 0.5} dB w.p. {0.2, 0.6, 0.2}, g⁽¹⁾ std 0.2 dB, ideal phase
correction. Unit tests: 100 trials per condition, noise 0–0.40 dB. E2: K = 64, γ ∈ {0.5, …, 0.99}. E5: γ = 0.9,
K ∈ {64, 128, 256, 512}. Both use the M2 simulator and a paper-faithful g⁽¹⁾ (< 0.1 Hz). Fig. 6: K = 256, 13 γ
values, dynamic types (i) and (ii). All runs: 500 trials, seed 42. Metrics: Eq. 29 χ (SNR as χ/(1−χ) and as the
paper's χ²/(1−χ²)), plus the RMS gain-tracking error in dB.

**Results and discussion.** Use §5.1–5.4 and §6. Key figures: `lambda_analysis.png` (cause),
`E5_plot.png` (fixed version now scales with K), `figure6_comparison.png` (paper vs ours, with the ideal-gain
ceiling).
