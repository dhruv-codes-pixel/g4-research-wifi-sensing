# M2 — Group 4 — Research Project: Wi-Fi sensing, CSI signal preprocessing

Base paper: V. V. Ratnam et al., "Optimal preprocessing of WiFi CSI for sensing applications," IEEE Trans. Wireless Commun., 23(9), 2024.

What M2 contains
- Report/  : M2 PDF report (and editable .docx)
- Code/    : CSI simulator, implementation of the paper's gain/phase correction + baselines, experiments E1-E5
- Results/ : figures E1-E5, results.json, summary_100_trials.txt
- Data/    : no external data at M2 (all data simulated, seeded)
- Video/   : link to walkthrough

How to run
    cd Code
    pip install -r requirements.txt
    cd src
    python run_experiments.py 100        # about 5 minutes; use 20 for a 1-minute test run
Outputs are written to ../../Results/
