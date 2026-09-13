# Milestone 1 — Proposal & Scoping

**Category:** Research in Wi-Fi Sensing
**Submission deadline:** Sep 7, 2026, 11:59 PM
**Status:** --

This folder holds everything submitted for M1: the proposal document, the individual reflections, and (once recorded) the milestone video.

## Contents

```
M1-proposal/
├── README.md                    
├── G4_PP2026_ECE310.pdf        
├── video/                       
└── reflections/                 
```

## Problem Statement

CSI-based Wi-Fi sensing methods have mostly been developed and validated on narrowband CSI under well-synchronized commodity devices. As sensing moves toward wider bandwidths and standards like **IEEE 802.11bf**, timing/synchronization errors and bandwidth variation are expected to affect signal fidelity — but how much, and at what point performance actually breaks, is not well characterized.

## Approach

Using the public **Widar 3.0** dataset, we:
1. Establish a baseline BVP classification accuracy on unperturbed CSI.
2. Inject controlled timing/synchronization offsets and sweep their magnitude.
3. Reduce CSI bandwidth/subcarrier count and sweep the reduction level.
4. Identify the knee/failure point on each sensitivity curve.
5. If degradation is significant, prototype a lightweight phase/timing correction and test recovery.

## Key Objectives

1. Measure how timing/synchronization offsets affect BVP classification accuracy.
2. Measure how reduced CSI bandwidth/subcarrier count affects BVP accuracy.
3. Identify the timing-error and bandwidth levels at which BVP performance degrades significantly (the knee/failure point).
4. If significant degradation is observed, develop and test a lightweight correction method to evaluate recoverability.

## Expected Outcomes

- Timing-error vs. BVP accuracy sensitivity curve
- Bandwidth vs. BVP accuracy sensitivity curve
- Knee/failure point identification for both axes
- A simple correction method (if feasible), with corrected vs. uncorrected comparison
- A reproducible experimental setup on public Widar 3.0 data, no new hardware
- A conclusion on whether/when sync error and bandwidth reduction break BVP generalization

## SOTA Anchoring

This work is positioned against recent Wi-Fi sensing generalizability surveys, cross-domain gesture recognition baselines (Widar3.0, SHARP), and clock-asynchronism/ISAC signal-processing literature, and is framed around the emerging **IEEE 802.11bf** sensing amendment — see [full references](../../docs/references.md).

## Dataset

**Widar 3.0** (IEEE DataPort, DOI [10.21227/7ZNF-QP86](https://ieee-dataport.org/open-access/widar-30-wifi-based-activity-recognition-dataset)) — CSI from Intel 5300 NICs across multiple environments, users, postures, orientations, and gesture classes. Used as the unperturbed baseline; perturbed variants are generated programmatically and documented in `data/perturbed/`.

## Deliverables

| Item | Location |
|------|----------|
| Proposal (PDF) | [`G4_PP2026_ECE310.pdf`](./G4_PP2026_ECE310.pdf) |

