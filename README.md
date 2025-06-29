# Muscle Synergy Analysis of Bipedal Squats under Varying Load Conditions

This repository contains MATLAB code to replicate the muscle synergy analysis conducted in the study:

> **"Investigating the impact of external load on muscle synergies during bipedal squats"**  
> Kenville et al., *European Journal of Applied Physiology*, 2024  
> DOI: [10.1007/s00421-024-05432-3](https://doi.org/10.1007/s00421-024-05432-3)

---

## 📈 Project Description

This project implements the full processing pipeline used to evaluate **how external loads affect muscle coordination patterns** (i.e., synergies) during the execution of **bipedal squats** using electromyography (EMG) data.

The study involved:
- 10 male participants
- Squats under 3 different loads (50%, 62.5%, 75% of 3-RM)
- Bilateral EMG recordings from 12 muscles

We replicate the synergy analysis using **Non-negative Matrix Factorization (NNMF)**, as described in the paper, and provide full preprocessing and plotting tools.

---

## 🧠 Methods

- **Preprocessing**
  - High-pass filter at 30 Hz
  - Hilbert rectification
  - Low-pass filter at 10 Hz
  - Normalization to max + standard deviation
  - Time normalization (resampled to 200 points)

- **Synergy Extraction**
  - NNMF using multiplicative update rules
  - Reconstruction accuracy threshold: 90%
  - Stop criteria: deltaRA < 3%
  - Final normalization of `W` and `H`

- **Plotting**
  - Temporal activation patterns (`H`) for all synergies across loads
  - Muscle contribution weights (`W`) with error bars and group-level bar graphs

---



