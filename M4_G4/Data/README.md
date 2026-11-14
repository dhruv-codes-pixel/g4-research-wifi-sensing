# Data

## Source

- Dataset: WiFi CSI Human Activities, folder `saved_csi_raw` (original CSI stage)
- URL: https://huggingface.co/datasets/aselya9185/wifi-csi-human-activities
- Original campaign: M. De Sanctis, R. Fallani, T. Rossi, E. Cianca, M. Ruggieri, V. Poulkov, "WiFi Sensing Through Digital Receive Beamforming and CSI," IEEE PIMRC 2025, DOI 10.1109/PIMRC62392.2025.11274992
- Hugging Face repository and processing code: A. Ussenova (MSc thesis, Tor Vergata, 2024/25), https://github.com/aselya9185/har-wifi-csi
- Licence: no formal licence is named. The dataset card states research and educational use and asks users to cite the PIMRC paper.
- Deviation from the project brief: the brief suggested IEEE DataPort or the SHARPax route. We used this open Hugging Face dataset because it is freely accessible without registration and contains labelled empty-room and activity recordings in two rooms.

## Contents

- 19 raw `.npy` recordings, 474 MB in total
- Data type `complex64`, shape `(4, N, 234)`; N (packets) varies per recording
- Two indoor offices (Room 1 `r1_*`, Room 2 `r2_*`)

| Class | Room 1 | Room 2 | Total |
|---|---|---|---|
| empty | 1 | 2 | 3 |
| sitting | 3 | 2 | 5 |
| standing | 3 | 2 | 5 |
| walking | 4 | 2 | 6 |

Full per-file details (class, session, shape, N, size, SHA-256, NaN/Inf count, zero-packet fraction, mean amplitude per antenna, estimated duration) are in `captures.csv`.

File naming in the source is inconsistent (`sit`/`sitting`, `walk`/`walking`, `r1_empty` without a number). `captures.csv` gives normalised classes and sessions; raw files keep their original names. `r1_walking_2` does not exist in the source dataset.

## Array layout

| Axis | Length | Meaning |
|---|---|---|
| 0 | 4 | receive antennas of the monitoring router |
| 1 | N | packets in time order |
| 2 | 234 | OFDM data subcarriers, ascending frequency |

Access pattern: `csi[antenna, packet, subcarrier]`.

## Tone map

- 80 MHz VHT channel, 256-point FFT, subcarrier spacing 312.5 kHz
- The 234 columns are FFT bins 6 to 250 without bins 25, 53, 89, 117, 127, 128, 129, 139, 167, 203, 231. Subcarrier index k = bin - 128, so k runs from -122 to +122 without the DC nulls (k = -1, 0, +1) and the pilots (k = +-11, +-39, +-75, +-103).
- Index per column: `docs/tone_map.csv` (column_index, subcarrier_index_k, fft_bin)
- Sub-band stitching boundary at k = -64 is column 56
- The upstream pipeline reconstructs the 8 pilots and k = +-1 by linear interpolation, giving 244 tones. That reconstruction is not part of this raw data.

## Capture details

- Hardware and firmware: Nexmon CSI on an Asus RT-AC86U router, IEEE 802.11ac VHT, 80 MHz, iPerf3 traffic between an access point and a client laptop
- Packet rate: not stored in the files. The thesis states nearly 5 packets per second (about ten minutes per activity). Recording lengths imply 7 to 17 minutes at 5 Hz, so treat 5 Hz as nominal and any time axis as approximate.
- Timestamps and sequence numbers: none
- Centre frequency and RF channel: not stated in the dataset card, the repository or the thesis
- Phase status: raw. The upstream pipeline applies PLL bias compensation and a linear regression transformation in later stages only. Evidence is in `docs/data_check.md`.
- Amplitude: uncalibrated units (about 1e12) with step-like level changes per antenna, consistent with AGC switching (not confirmed by the source)

## Raw data policy

- Raw files are never modified or renamed.
- Raw files are not committed to the repository because of their size.
- Integrity is verified with the SHA-256 checksums in `docs/raw_sha256.txt`.
- All derived outputs go to `Results/`, never into the raw folder.

## How to obtain the data

```
pip install huggingface_hub numpy matplotlib
```
```
python -c "from huggingface_hub import snapshot_download; snapshot_download('aselya9185/wifi-csi-human-activities', repo_type='dataset', allow_patterns='saved_csi_raw/*', local_dir='.')"
```
Verify on Windows PowerShell:
```
Get-FileHash saved_csi_raw\*.npy -Algorithm SHA256
```
Compare the hashes with `docs/raw_sha256.txt`.

## Loading and inspection code

- `csi_loader.py`: `load_csi(path)` returns a capture with the `(4, N, 234)` array, room, class, session and the k index; `layout` can reorder axes for code that expects a different layout
- `inspect_npy.py`: generates `captures.csv`, `docs/data_check.md`, `docs/tone_map.csv`, `docs/raw_sha256.txt` and the diagnostic plot

```
python inspect_npy.py --data saved_csi_raw --out docs --csv captures.csv
```

## Known limitations

- No timestamps; durations are estimates from a nominal 5 Hz rate
- Amplitude gain steps (likely AGC) must be corrected before cross-packet comparison
- Centre frequency not documented
- Few recordings per class (3 to 6), and `r1_walking_2` is absent in the source
- The raw-phase verdict was measured on `r1_empty` and checked against the thesis description, not verified per file

## Related files

- `docs/data_check.md`: full data check with axes, tone map, raw vs sanitised phase, amplitude, timestamps
- `captures.csv`: per-recording summary
- `docs/tone_map.csv`: subcarrier index per column
- `docs/raw_sha256.txt`: checksums of the raw files
