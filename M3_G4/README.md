# M3: Wi-Fi CSI Preprocessing Under Realistic TGn Channels, Ultra-Wide Bandwidth, and Hardware Imperfections
**Group G4 — ECE 310: Wireless Communications Major Project (Milestone 3)**

[![Milestone](https://img.shields.io/badge/Milestone-M3%20Evaluation-blue.svg)](Report/M3_g4-research-wifi-sensing_Report.pdf)
[![Python](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11-brightgreen.svg)](Code/)
[![Status](https://img.shields.io/badge/Status-Complete-success.svg)](Protocol.md)
[![Target](https://img.shields.io/badge/Target-IEEE%20Transactions-orange.svg)](Report/M3_g4-research-wifi-sensing_Report.pdf)

---

### Executive Overview

This repository houses the final Milestone 3 (M3) code, experimental results, and technical documentation for Group G4. Our research investigates the robustness of optimal Wi-Fi Channel State Information (CSI) gain and phase preprocessing algorithms, building upon the theoretical foundation of **Ratnam et al. (IEEE Transactions on Wireless Communications, Sept. 2024)**.

Milestone 3 extends the preliminary M2 findings by:
1. **Upgrading Channel Realism:** Transitioning from an idealized 6-tap synthetic channel to the standardized **IEEE 802.11 TGn Model C tapped-delay-line (TDL)** channel (14 taps, 33.44 ns RMS delay spread, dual overlapping Rayleigh fading clusters).
2. **Exhaustive Stress-Testing (500 Monte Carlo Trials per Point):**
   - **E1 (Static Power Fraction):** Evaluating phase recovery across $\gamma \in [0.5, 0.99]$.
   - **E3 (Synchronization Stress):** Stressing timing offsets up to $\tau_{\text{max}} = 3.2\;\mu\text{s}$ ($4\times$ the physical cyclic prefix length).
   - **E4 (Bandwidth Scaling):** Scaling subcarrier counts across $K \in \{64, 128, 256, 512\}$, modeling channels from 20 MHz to 160 MHz in Wi-Fi 6/6E/7.
3. **Downstream Sensing Application Validation:** Injecting a 20 Hz Doppler velocity target and evaluating peak estimation accuracy and line-to-floor ratios with paired Wilcoxon signed-rank significance testing and bootstrap confidence intervals.
4. **Real-World CSI Pipeline Engineering:** Constructing an end-to-end ingestion and preprocessing engine (`Code/real_data.py`) for multi-device commercial hardware captures, accompanied by a documented audit of remote dataset hosting access barriers (4TU and Zenodo).

---


---

### Team Members & Contributions

| Member | Role | Primary Responsibilities |
|---|---|---|
| **Pratyusha** | Team Lead & Systems Architect | Evaluation protocol, journal positioning, report integration, repository management |
| **Jay Patel** | Realistic Channel Modeling | TGn Model C implementation (`channel_tgn.py`), E1/E3/E4 batch experiments, Fig 7 overlay |
| **Dhruv** | Hardware CSI & Dataset Audit | Real-data pipeline (`real_data.py`), 4TU/Zenodo access audit (`data_check.md`, `log.txt`) |
| **Lay** | Downstream Sensing & Statistics | Doppler sensing experiment (`doppler_exp.py`), rSVD extraction, Wilcoxon tests (`stats_utils.py`) |
| **Member B** | Gain Dynamics & Complexity | Algorithm 3 anomaly diagnosis, discrete AGC evaluation, algorithmic profiling across $K$ |

For full operational procedures, task allocation rationale, and our long-term roadmap toward journal submission, see **[Protocol.md](Protocol.md)**.

---

### Key Technical Findings

1. **Realistic Multipath Exposes Baseline Vulnerability:**
   Under the standardized TGn Model C channel, the conventional least-squares phase unwrap baseline (Eq. 21) collapses from $+15.2\text{ dB}$ (in M2) down to **$-1.52\text{ dB}$** at $\gamma = 0.9$. In contrast, the proposed Algorithm 4 + Eq. (25) remains completely impervious, achieving **$16.53\text{ dB}$** (virtually identical to M2's $16.50\text{ dB}$).
2. **Bandwidth Scaling Strongly Favors Algorithm 4:**
   As channel bandwidth expands from $K = 64$ (20 MHz) to $K = 512$ (160 MHz), the proposed method **gains $+4.55\text{ dB}$** (reaching $21.27\text{ dB}$), while Eq. (21) and Eq. (22) degrade by $-23.86\text{ dB}$ and $-10.83\text{ dB}$, respectively. At $K = 512$, the proposed method outperforms Eq. (21) by **$+45.05\text{ dB}$**.
3. **Synchronization Invariance:**
   The proposed method demonstrates flat performance ($16.42\text{ dB}$ to $16.58\text{ dB}$) across timing errors ranging from 100 ns to 3200 ns, showing extreme resilience even when timing errors exceed the 0.8 $\mu\text{s}$ cyclic prefix.
4. **Direct Translation to Sensing Fidelity:**
   In 20 Hz Doppler velocity sensing, Algorithm 4 consistently yields the sharpest spectral peaks and lowest frequency estimation error, with paired Wilcoxon signed-rank tests confirming statistically significant superiority ($p < 0.001$) over uncorrected and baseline methods.

---


