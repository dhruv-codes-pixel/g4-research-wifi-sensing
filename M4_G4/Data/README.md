# 20 MHz 802.11ax Wi-Fi CSI Data (split archives)

This repository holds split archives of the **20 MHz** part of a public Wi-Fi CSI dataset. The original `20MHz.zip` (~481 MB) is too large for a single GitHub upload, so it was divided into the smaller `ZIP*.zip` files here.

## Original dataset

- **Title:** *Towards a Quantitative Analysis of CSI for AI/ML Based Sensing* (v2: *802.11ax Wi-Fi CSI Measurements for Wireless Sensing*)
- **Authors:** Elena Tonini, Renato Lo Cigno (University of Brescia); v1 also Francesco Gringoli, Marco Cominelli
- **Source:** Zenodo – https://zenodo.org/records/14651521
- **DOI:** [10.5281/zenodo.14651521](https://doi.org/10.5281/zenodo.14651521)
- **Licence:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
- **Paper:** E. Tonini et al., "Towards a Quantitative Analysis of CSI for AI/ML Based Sensing," IEEE WCNC 2025, doi: [10.1109/WCNC61545.2025.10978183](https://doi.org/10.1109/WCNC61545.2025.10978183)

Please cite the paper and the Zenodo DOI (check the record page for the exact version DOI of your copy). The data is redistributed unchanged under the original licence; all credit belongs to the authors.

## Files NOT attached

Two archives could not be attached because they exceed the file-size limit for upload:

| File | Content |
|---|---|
| `ZIP__20mhz__2_.zip` (~88 MB) | `capture0–3.pcap` (empty room, 2024-07-12) |
| `ZIP_mat.zip` (~85 MB) | `capture0–3.mat` (empty room, 2024-07-12) |

Download them (and the rest of the dataset) from the Zenodo link above.

## What is attached

The small remaining parts contain only a `config.json` and empty macOS `._*` metadata stubs (safe to ignore). The folder structure of the original archive was **not** preserved.

## Using the data

1. Download the full 20 MHz data from Zenodo and unzip it, keeping this layout:
   ```
   20MHz/0people/<date>/  1people/<date>/  4people/<date>/   (captureN.mat, captureN.pcap, config.json)
   ```
2. `.mat` files are MATLAB v7.3 (HDF5): variable `capture`, 256 FFT bins, 242 active tones (−122…−2, +2…+122).
3. Point the code to it:
   ```bash
   python run_m4_experiment.py --data-root /path/to/20MHz --captures captures.csv --out results/real_20mhz
   ```
