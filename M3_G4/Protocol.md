# ECE 310 — Wireless Communications: Major Project Milestone 3 (M3)
**Project Title:** How Well Does Wi-Fi CSI Gain and Phase Correction Hold Up as Bandwidth and Synchronization Error Change?  
**Base Reference:** V. V. Ratnam, H. Chen, H.-H. Chang, A. Sehgal, and J. Zhang, *"Optimal Preprocessing of WiFi CSI for Sensing Applications,"* IEEE Transactions on Wireless Communications, vol. 23, no. 9, pp. 10820–10833, Sept. 2024.  
**Target Publication Venue:** IEEE Transactions on Wireless Communications / IEEE Transactions on Mobile Computing (Regular Research Article)  
**Lead Author / Integration Lead:** Pratyusha 

---

### Executive Summary & Research Vision

Milestone 3 (M3) marks the formal transition from initial technical implementation (M2) to exhaustive experimental stress-testing, downstream sensing validation, and empirical State-of-the-Art (SOTA) positioning. While M2 established baseline reproducibility over an idealized 6-tap Rayleigh fading model, the unsaid strategic goal of our 5-member team is to formulate an archival-quality study suitable for submission to a flagship journal (**IEEE Transactions on Wireless Communications** or **IEEE Transactions on Mobile Computing**). 

Our central scientific thesis addresses a critical unanswered question in RF sensing:
> *"Does optimal Wi-Fi CSI gain and phase preprocessing maintain its theoretical superiority when subjected to standardized frequency-selective multipath (IEEE 802.11 TGn Model C), extreme timing synchronization errors ($\tau \gg T_{\text{CP}}$), and modern ultra-wide bandwidth scaling ($K = 64 \to 512$ subcarriers, representing 20 MHz to 160 MHz channels in Wi-Fi 6/6E/7)?"*

To address this question rigorously, the work was structured across five specialized roles, ensuring complete accountability, incremental Git commit tracking, and verifiable experimental artifacts.

---

#### Member 1: Pratyusha (Team Lead, Framework Synthesis & Journal Positioning)
* **Assigned Scope:** System-level architecture,evaluation protocol definition, integration of all members' findings into the master manuscript, and journal-level positioning.
* **Key Tasks Executed:**
  1. **Evaluation Protocol Definition:** Established the standardized testing matrix (500 Monte Carlo trials per parameter point, 95% bootstrap confidence intervals for median estimates, and non-parametric paired significance tests).
  2. **Real-World Evidence Management:** Analyzed the root-cause connectivity barriers encountered when fetching external datasets (4TU ResearchData server HTTP 403 refusal and Zenodo API blocks); formulated the transparent disclosure policy confirming zero data falsification.
  3. **Journal-Level Manuscript Composition:** Authored the comprehensive M3 Report following IEEE Transactions guidelines, synthesizing channel validation, phase stress tests, bandwidth expansion dynamics, and downstream Doppler performance.
* **Proof of Work / Required Deliverables in Submission:**
  - Complete, unified `M3_g4-research-wifi-sensing_Report.pdf` with publication-grade figures, detailed mathematical proofs, and critical analyses.
  - Master `Protocol.md` detailing operational division and workflow execution.
  - Overall repository cleanliness, module modularity, and reproducible execution entry points.

---

#### Member 2: Jay Patel (Realistic Channel Modeling & Core Phase Experiments)
* **Assigned Scope:** Replacing the synthetic 6-tap exponential channel of M2 with the standardized IEEE 802.11 TGn Model C tapped-delay-line (TDL) profile; scaling Monte Carlo trials from 100 to 500; evaluating phase-correction baselines vs. the proposed method.
* **Key Tasks Executed:**
  1. **Standardized Channel Implementation (`channel_tgn.py`):** Encoded the 14-tap Model C delays (0 to 200 ns) and dual overlapping Rayleigh cluster profiles (Cluster 1: taps 1–10; Cluster 2: taps 7–14), combining overlapping tap powers in the linear domain before complex tap drawing.
  2. **Zero-Slope Phase Preservation (`csi_sim.py`):** Enforced that the frequency-selective realization $b_k$ maintains $\sum \angle(b_k b_{k+1}^*) = 0$ and $\frac{1}{K}\sum |b_k|^2 = \gamma$, preventing synthetic linear slopes from biasing timing estimators.
  3. **Batch Monte Carlo Runner (`run_m3.py`, `run_single_curve.py`):** Executed 500 trials per point across:
     - **E1:** $\gamma \in [0.5, 0.99]$ at $K = 64$ (phase recovery vs. static power fraction).
     - **E3:** $\tau_{\text{max}} \in [100\text{ ns}, 3200\text{ ns}]$ at $K = 64, \gamma = 0.9$ (timing synchronization stress).
     - **E4:** $K \in \{64, 128, 256, 512\}$ subcarriers (bandwidth scaling from 20 MHz to 160 MHz).
  4. **Empirical Overlay with Published Art:** Digitized Ratnam et al. Fig. 7(a) markers to generate an empirical comparison plot between TGn Model C simulation and published curves.
* **Proof of Work / Required Deliverables in Submission:**
  - Source modules: `channel_tgn.py`, `csi_sim.py`, `run_m3.py`, `run_single_curve.py`.
  - Raw JSON execution data: `M3_results.json` containing median SNR, Q1/Q3, and 95% bootstrap CIs.
  - Result Figures: `M3_channel_PDP.png`, `M3_channel_frequency_response.png`, `M3_old_vs_new_channel.png`, `M3_E1_phase_vs_gamma.png`, `M3_E3_phase_vs_timing.png`, `M3_E3_CP_ISI_diagnostic.png`, `M3_E4_phase_vs_K.png`, `M3_Figure7_overlay.png`.

---

#### Member 3: Dhruv (Real-World Dataset Acquisition & Hardware-In-The-Loop Validation)
* **Assigned Scope:** Bridging synthetic simulations with real-world Wi-Fi CSI hardware captures; procuring raw multi-device CSI traces from published repositories; constructing an end-to-end preprocessing validation pipeline.
* **Key Tasks Executed:**
  1. **Dataset Identification & Specification (`data_check.md`, `dataset_info.csv`):**
     - Target 1: Portner et al. (Sensession, TU Delft / 4TU.ResearchData, DOI: `10.4121/48e1331b-91df-4e1a-906b-b0a0e9d16c74`) comprising 8 receivers (USRP N2954-R, Intel AX210, Intel IWL5300, Qualcomm QCA AR9462, Asus RT-AC86U Nexmon, ESP32-S3) connected via coaxial cables and equal-length power splitters. Target file: `base_campaign_ch11.zip` (539 MiB, 10,000 packets).
     - Target 2 (Fallback): OpenCSI (Zenodo Record `19861032`, LCN 2026) featuring 8 ESP32 nodes over 802.11n/ax links.
  2. **Access Barrier Documentation (`log.txt`):** Systematically documented every server-level connection refusal encountered during fetch operations (HTTP 403 tunnel forbidden errors on 4TU, API handshake failures, and DOI browser timeouts).
  3. **Hardware Preprocessing Engine (`real_data.py`):**
     - Built an automated pipeline capable of ingesting nested parquet arrays (`csi.parquet` + `meta.parquet`).
     - Handled device-specific quirks: AGC warm-up packet truncation (first 200 frames), null subcarrier excision, DC removal, and non-uniform subcarrier spacing safeguards (automatically skipping Eq. 22 / Algo 4 on 30-subcarrier IWL5300 captures).
     - Formulated amplitude and phase residual metrics relative to the static time-averaged channel.
     - Validated pipeline functionality via synthetic traces in native parquet schema (pipeline test verified, output discarded).
* **Proof of Work / Required Deliverables in Submission:**
  - Automated pipeline: `real_data.py`.
  - Comprehensive access audit: `data_check.md` and `Results/real_data/log.txt`.
  - Dataset registry: `Results/real_data/dataset_info.csv`.
  - Status manifest: `Results/real_data/experiment_config.json` (explicitly recording `status: not_run` with complete parameters).

---

#### Member 4: Lay (Downstream Sensing Applications & Non-Parametric Statistical Inference)
* **Assigned Scope:** Proving that CSI preprocessing directly translates to tangible gains in downstream sensing applications (specifically Doppler velocity estimation); establishing rigorous statistical significance for all observed improvements.
* **Key Tasks Executed:**
  1. **Doppler Sensing Simulation Engine (`doppler_exp.py`):**
     - Synthesized multi-frame CSI ($P = 512$ frames, $T_{\text{rep}} = 5\text{ ms}$) impinging on a moving reflector inducing a 20 Hz Doppler frequency shift.
     - Implemented randomized Singular Value Decomposition (rSVD) to isolate dominant dynamic multipath reflections without requiring ground-truth reflector channel knowledge.
     - Tested 5 correction regimes: No Correction, Eq. (6) + Eq. (22), Eq. (6) + Eq. (21), Eq. (6) + Algorithm 4, and Ideal Oracle Correction.
     - Evaluated across 6 configurations: $K \in \{64, 256, 512\} \times \tau_{\text{max}} \in \{100\text{ ns}, 800\text{ ns}\}$, with 200 independent trials per configuration ($6 \times 200 \times 5$ tensor).
  2. **Statistical Inference Utilities (`stats_utils.py`):**
     - Executed 10,000 percentile-bootstrap resamples to establish non-parametric 95% confidence intervals for median estimates.
     - Conducted paired two-sided Wilcoxon signed-rank tests comparing each correction strategy against the uncorrected baseline.
     - Applied Holm-Bonferroni step-down adjustments to control family-wise error rates ($\alpha = 0.05$).
  3. **Caveat & Estimation Boundary Analysis:** Documented the FFT interpolation limit, noting that zero-padded spectral peak estimation quantizes to 0.0195 Hz relative to the native 0.39 Hz bin spacing.
* **Proof of Work / Required Deliverables in Submission:**
  - Source modules: `doppler_exp.py`, `stats_utils.py`.
  - Processed data arrays & summaries: `results/doppler_results.npz`, `doppler_summary.csv`, `doppler_bootstrap_summary_200trials.csv`, `doppler_wilcoxon_tests_200trials.csv`.
  - Spectral plots: Doppler peak error, line-to-floor ratio, and 6 spectrum figures across $K$ and $\tau$.

---

#### Member 5: Member B (Gain Dynamics, AGC Profiling & Algorithmic Complexity)
* **Assigned Scope:** Profiling the gain-correction pipeline; analyzing the interaction between discrete AGC switching and frequency-selective fading; isolating failure modes in the maximum-likelihood AGC estimator (Algorithm 3).
* **Key Tasks Executed:**
  1. **Gain Pipeline Benchmarking:** Evaluated Eq. (6) (RMS normalization), Algorithm 1 (DBSCAN clustering of power increments with 0.1 Hz moving average filter), and Algorithm 3 (ML line search over candidate AGC step sizes $\lambda$).
  2. **Algorithm 3 Anomaly Diagnosis:** Identified and documented a numerical trap in the M2 implementation of Algorithm 3, where the objective function $\sigma^2 + \lambda^2 \mathcal{D}(\lambda/\sigma)$ converges to $\lambda \approx 3.5\text{ dB}$ instead of the true 0.5 dB step under TGn multipath, causing residual amplitude variance.
  3. **Computational & Memory Complexity Analysis:** Formulated the asymptotic complexity profiles across subcarrier count $K$:
     - Baseline Eq. (21): $\mathcal{O}(P \cdot K)$ least-squares unwrap (vulnerable to phase discontinuities).
     - Baseline Eq. (22): $\mathcal{O}(P \cdot K)$ adjacent-subcarrier correlation (limited resolution).
     - Algorithm 4 + Eq. (25): $\mathcal{O}(P \cdot K)$ vectorized 7-tap sliding window robust unwrap + weighted least-squares closed form.
  4. **Hardware Footprint Estimation:** Computed memory bounds showing that vectorized Algorithm 4 requires under 1.8 MB of RAM at $K = 512, P = 300$, rendering it fully deployable on embedded edge platforms (e.g., ESP32-S3 or OpenWrt routers).
* **Proof of Work / Required Deliverables in Submission:**
  - Gain evaluation scripts and mathematical diagnosis of Algorithm 3 grid convergence.
  - Tabulated computational complexity and memory footprint comparison across $K \in \{64, 128, 256, 512, 1024\}$.
  - Technical analysis integrated directly into Section 2 and Section 8 of the Master Report.

---

### 2. Strategic Research Roadmap: Expanding to Ultra-Wide Bandwidth & Journal Publication

To fulfill our team's vision of publishing these findings in **IEEE Transactions on Wireless Communications**, our post-M3 research trajectory focuses on three critical pillars:

```
+----------------------------------------------------------------------------------------------------+
|                         POST-M3 EXPANSION ROADMAP TOWARD IEEE TRANSACTIONS                         |
+----------------------------------------------------------------------------------------------------+
|  Phase 1: Ultra-Wide Bandwidth Scaling (Wi-Fi 7 / 320 MHz)                                         |
|  - Extend subcarrier evaluation from K = 512 (160 MHz) to K = 1024 and K = 2048 (320 MHz).         |
|  - Model punctured subcarrier patterns (OFDMA RU allocations in 802.11ax/be).                      |
|  - Quantify the asymptotic SNR saturation limit of Algorithm 4 as K -> infinity.                   |
+----------------------------------------------------------------------------------------------------+
|  Phase 2: Over-The-Air Hardware Testbed with Dual USRP X310                                        |
|  - Overcome public repository hosting barriers by capturing native CSI via NI USRP X310 SDRs.      |
|  - Validate coaxial baseline (zero-multipath ground truth) vs. indoor rich-multipath environment.   |
|  - Compare hardware-measured dynamic CSI SNR against our TGn Model C simulation curves.            |
+----------------------------------------------------------------------------------------------------+
|  Phase 3: Multi-Antenna MIMO & Spatio-Temporal Preprocessing                                       |
|  - Generalize Algorithm 4 to 2x2 and 4x4 MIMO cross-antenna phase center variations.               |
|  - Derive the Cramér-Rao Lower Bound (CRLB) for joint gain-timing-phase estimation under TGn.      |
+----------------------------------------------------------------------------------------------------+
```

---

### 3. Verification Checklist for Milestone Evaluation

| Rubric Evaluation Item | Status | Primary Evidence in Repository | Responsible Member |
|---|---|---|---|
| **Final Proposed Method Established** | **COMPLETE** | Fully vectorized Algorithm 4 + Eq. (25) integrated with TGn Model C in `Code/csi_sim.py` | Pratyusha & Jay Patel |
| **Experimental Setup Documented** | **COMPLETE** | Full specification of TGn Model C PDP, noise parameters, $\tau$ bounds, and bootstrap protocols | Jay Patel & Pratyusha |
| **Demonstrated Effectiveness** | **COMPLETE** | E1, E3, E4 results exhibiting 16.53 dB median SNR at $\gamma=0.9$, $+4.5\text{ dB}$ gain across $K$ | Jay Patel |
| **Baseline Comparisons** | **COMPLETE** | Comprehensive comparison vs. Eq. (21) and Eq. (22), showing Eq. (21) collapses by up to 34 dB | Jay Patel & Member B |
| **Downstream Sensing Proof** | **COMPLETE** | Doppler velocity experiment confirming 20 Hz peak recovery with minimal error and statistical testing | Lay |
| **Real-Data Pipeline & Audit** | **COMPLETE** | `real_data.py` operational; transparent audit of 4TU and Zenodo access barriers in `log.txt` | Dhruv |
| **Revised SOTA Positioning** | **COMPLETE** | Detailed SOTA revision proving Algorithm 4's bandwidth scaling superiority for Wi-Fi 6E/7 | Pratyusha |
| **Submission Structure Adherence** | **COMPLETE** | Formatted strictly as `M3_g4-research-wifi-sensing/` with `Report/` containing publication PDF | Pratyusha |


