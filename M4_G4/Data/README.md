# 20 MHz 802.11ax Wi-Fi CSI Data (split archives)

This repository holds split archives of the **20 MHz** part of a public Wi-Fi CSI dataset. The original `20MHz.zip` (~481 MB; 504.8 MB on Zenodo) is too large for a single GitHub upload, so it was divided into the smaller `ZIP*.zip` files here.

## Original dataset

- **Title:** *802.11ax Wi-Fi Channel State Information (CSI) Measurements for Wireless Sensing* (Zenodo, Version 2, published October 2025; Version 1 was titled *Towards a Quantitative Analysis of CSI for AI/ML Based Sensing*)
- **Authors:** Elena Tonini, Renato Lo Cigno, Francesco Gringoli (University of Brescia), Marco Cominelli
- **Source:** Zenodo - https://zenodo.org/records/17225694
- **DOI (this version, v2):** [10.5281/zenodo.17225694](https://doi.org/10.5281/zenodo.17225694)  |  all-versions DOI: [10.5281/zenodo.14651521](https://doi.org/10.5281/zenodo.14651521)
- **Licence:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
- **Associated paper (v1):** E. Tonini, F. Gringoli, R. Lo Cigno, M. Cominelli, "Towards a Quantitative Analysis of CSI for AI/ML Based Sensing," IEEE WCNC 2025, doi: [10.1109/WCNC61545.2025.10978183](https://doi.org/10.1109/WCNC61545.2025.10978183)
- **Associated paper (v2):** E. Tonini, R. Lo Cigno, "Leveraging Mutual Information in Stochastic CSI Analysis for Wi-Fi Sensing," IEEE/IFIP WONS 2026, pp. 97-104

**Verified:** the 20 MHz archive used here has the same size (504,837,103 bytes) and MD5 (`311b862e49952929c162092e469dc691`) as `20MHz.zip` on the Zenodo v2 record; the dataset README matches too (MD5 `9052511dcf4490d12599f3ba16b78b0d`).

Please cite the dataset DOI and the associated paper(s). The data is redistributed unchanged under the original licence; all credit belongs to the authors.

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
