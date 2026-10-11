# M2 Code — Reproduction of Ratnam et al. (Samsung), IEEE TWC 2024 (simulation arm)

Files
- `src/csi_sim.py`        : CSI simulator (paper eq. 1-2) + gain correction (Eq.6, Algo 1, Algo 3) + phase correction (Eq.21, Eq.22, Algo 4 + Eq.25) + metric (Eq.29)
- `src/run_experiments.py`: runs experiments E1-E5, saves plots + results.json to `../Results/`

Run
    pip install -r requirements.txt
    cd src && python run_experiments.py 100      # 100 trials per point (~5 min); paper uses 2000

Known simplifications vs. paper (state these in the report)
1. Static channel b_k = random exponential-PDP multipath (NOT the 802.11ax Model-C channel).
2. Only dynamic-component type (i) (iid complex Gaussian) is simulated; type (ii) is not.
3. 100 trials/point instead of 2000.
4. K=64, subcarrier spacing 312.5 kHz (Ts=3.2 us). Paper also uses K=256 with the same Ts.
5. Robust unwrap Eq.(26) and Eq.(16) implemented from the arXiv text; Eq.(16) rounding form should be re-checked against the published PDF.
6. Sign convention in Eq.(21) fixed to be consistent with Eq.(2) (fit -angle(h) = 2*pi*f*tau + psi).
7. Metric: chi (Eq.29) is a squared correlation; SNR = chi/(1-chi). The paper's exact SNR scale is not stated, so compare TRENDS/ordering, not absolute dB.

Sanity check (compares estimators with the true simulated values)
    cd src && python sanity_checks.py
Known open issue: our Algorithm 3 picks lambda ~ 3.5 dB instead of the true 0.5 dB AGC step, so it does not track the gain (see report Sec. 6.4). To be fixed for M3.
