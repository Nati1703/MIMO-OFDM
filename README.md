# MIMO-OFDM Channel Estimation and Detection: Classical and Learned Approaches

## Overview

End-to-end MIMO-OFDM link-level simulator for evaluating channel estimation and signal detection algorithms under realistic fading conditions. Compares classical approaches (LS, MMSE, DFT-based estimation; ZF, MMSE, OSIC detection) with a neural network-based channel estimator.

This project integrates concepts from coursework in Statistical Signal Processing, Adaptive & Array Signal Processing, and Numerical Linear Algebra for Signal Processing at the Technical University of Munich (TUM).

## System Architecture

```
TRANSMITTER
────────────────────────────────────────────────────────
Bits → QAM Mapper → S/P
                                            ↓
                            Antenna 1: IFFT → Add CP → Transmit
                            Antenna 2: IFFT → Add CP → Transmit
                                  ⋮

CHANNEL
────────────────────────────────────────────────────────
  Frequency-selective MIMO channel (Kronecker correlation model)
  H[k] = per-subcarrier Nr × Nt channel matrix
  y[k] = H[k] · x[k] + n[k]  at each subcarrier k

RECEIVER
────────────────────────────────────────────────────────
  Rx Ant 1: Remove CP → FFT → Y₁[k]  ─┐
  Rx Ant 2: Remove CP → FFT → Y₂[k]  ─┤
       ⋮                                │
                                        ↓
  ┌─────────────────────────────────────────────┐
  │  Per subcarrier k = 0, 1, ..., N-1:         │
  │                                             │
  │  1. Channel Estimation                      │
  │     Ĥ[k] from pilots using LS / MMSE /     │
  │     DFT-based / Neural Network              │
  │                                             │
  │  2. MIMO Detection                          │
  │     x̂[k] from y[k] and Ĥ[k] using        │
  │     ZF / MMSE / OSIC                        │
  └─────────────────────────────────────────────┘
                        ↓
  P/S → QAM Demap → Deinterleaver → FEC Decode → Bits
```

## Modules

### Module 1: MIMO-OFDM Baseline
- **Channel generation:** i.i.d. Rayleigh and spatially correlated (Kronecker model) MIMO channels
- **Frequency selectivity:** exponential power delay profile, configurable delay spread, FFT across taps for per-subcarrier H[k]
- **OFDM transceiver:** IFFT modulation, cyclic prefix insertion/removal, FFT demodulation
- **Channel estimation:** Least Squares (LS) and Minimum Mean Square Error (MMSE)
- **MIMO detection:** Zero-Forcing (ZF) and MMSE
- **Verification:** BER curves vs analytical performance in flat Rayleigh fading

### Module 2: Enhanced Algorithms
- **DFT-based channel estimation:** IFFT → retain L significant taps → zero noise-only region → FFT. Reduces noise by factor L/N without needing channel statistics.
- **OSIC detection:** Ordered Successive Interference Cancellation with SINR-based ordering. Detects strongest stream first, subtracts, repeats.
- **BER comparison:** LS vs MMSE vs DFT-enhanced LS; ZF vs MMSE vs OSIC across SNR range.

### Module 3: Learned Channel Estimator
- **Architecture:** 3-layer fully connected neural network (MLP) in PyTorch
- **Input:** received pilot signals Y_p (complex-valued, split into real/imaginary)
- **Output:** channel estimate Ĥ at all subcarriers
- **Training data:** generated from Module 1 MATLAB simulator (channel realizations + noisy observations)
- **Key result:** achieves near-MMSE MSE performance without requiring the channel correlation matrix R_HH or noise variance — the network learns these implicitly from data

## Project Structure
```
mimo-ofdm/
│
├── README.md
│
├── matlab/
│   ├── channel/
│   │   ├── Rayleigh_model.m       % Base SISO Rayleigh generator
│   │   ├── Rician_model.m         % Base SISO Rician generator
│   │   ├── rayleigh_mimo.m        % i.i.d. Rayleigh MIMO channel
│   │   ├── rician_mimo.m          % Rician MIMO with configurable K-factor
│   │   ├── kronecker_mimo.m       % Spatially correlated MIMO (Kronecker)
│   │   ├── exp_pdp.m              % Exponential power delay profile
│   │   ├── freq_selective_mimo.m  % Static frequency-selective MIMO channel
│   │   ├── doppler_spectrum.m     % Clarke/Jakes Doppler spectrum
│   │   ├── fwgn_model.m           % Time-correlated fading (FWGN)
│   │   └── time_varying_mimo.m    % Doubly-dispersive MIMO channel
│   │
│   ├── ofdm/
│   │   ├── ofdm_mod.m            % IFFT + CP insertion
│   │   ├── ofdm_demod.m          % CP removal + FFT
│   │   └── qam_mod_demod.m       % QAM mapping/demapping
│   │
│   ├── estimation/
│   │   ├── ls_ce.m               % LS channel estimation
│   │   ├── mmse_ce.m             % MMSE channel estimation
│   │   └── dft_ce.m              % DFT-based denoised estimation
│   │
│   ├── detection/
│   │   ├── zf_detector.m         % Zero-forcing detection
│   │   ├── mmse_detector.m       % MMSE detection
│   │   └── osic_detector.m       % OSIC with SINR-based ordering
│   │
│   ├── preliminary/
│   │   ├── dft_review.m          % DFT/DSP fundamentals verification
│   │   └── channel_verification.m % Channel models verification
│   │
│   ├── simulate.m                % Main simulation script (BER curves)
│   └── gen_training_data.m       % Generate (Y_pilot, H) pairs for NN training
│
├── python/
│   ├── learned_estimator.py      % PyTorch neural network estimator
│   ├── train.py                  % Training loop
│   ├── evaluate.py               % Load results, compare with classical, plot
│   └── data/                     % Saved .mat training/test data
│
├── results/
│   ├── ber_estimation.png        % BER: LS vs MMSE vs DFT vs NN
│   ├── ber_detection.png         % BER: ZF vs MMSE vs OSIC
│   └── mse_comparison.png        % MSE vs SNR for all estimators

```
## Channel Models

### Flat Rayleigh Fading
```
H = (1/√2)(randn(Nr,Nt) + j·randn(Nr,Nt))
```
Each entry is CN(0,1). Same H for all subcarriers.

### Kronecker Correlation Model
```
H = R_rx^(1/2) · H_iid · R_tx^(1/2)
```
Introduces spatial correlation from antenna spacing and angular spread.

### Frequency-Selective Fading
```
h[l] ~ CN(0, σ²_l)    for l = 0, ..., L-1
σ²_l = exp(-l/τ_rms)   (exponential PDP)
H[k] = FFT(h, N)       (per-subcarrier channel)
```

## Key Algorithms

### Channel Estimation

| Method | Formula | Needs | MSE |
|--------|---------|-------|-----|
| LS | Ĥ = Y/X (at pilots) | Nothing | σ²/P_pilot |
| MMSE | Ĥ = R_HH(R_HH + (1/SNR)·I)⁻¹·Ĥ_LS | R_HH, σ² | Lower than LS |
| DFT-based | IFFT → keep L taps → FFT | L only | ≈ LS · (L/N) |
| Neural Net | Ĥ = f_θ(Y_pilot) | Training data | ≈ MMSE |

### MIMO Detection

| Method | Performance    |
|--------|----------------|
| ZF     | Baseline       |
| MMSE   | Better than ZF |
| OSIC   | Successive detection + cancellation | & ML Approaches  |


Planned Extensions: Iterative Receivers & Modern Channel Coding
To further bridge classical communication theory with modern learning-based approaches, the repository will be expanded to include iterative receiver architectures informed by concepts from the Modern Channel Codes lecture at Technical University of Munich.
Motivation
Current receiver pipeline:
Channel Estimation → Detection → Decoding   (one-shot)
Planned extension:
Channel Estimation ⇄ Detection ⇄ Decoding   (iterative / turbo processing)


### Python (Module 3)
```bash
cd python/
pip install torch scipy matplotlib

# Train the neural network estimator
python train.py

# Evaluate and compare with classical methods
python evaluate.py    # produces comparison plots in results/
```

## Dependencies

- **MATLAB** R2020a or later
- **Python** 3.8+
- **PyTorch** 1.9+
- **SciPy** (for loading .mat files)
- **Matplotlib** (for plotting)

## References

1. Y.S. Cho, J. Kim, W.Y. Yang, C.G. Kang, *MIMO-OFDM Wireless Communications with MATLAB*, Wiley, 2010.
2. Tzi-Dar Chiueh, Pei-Yun Tsai Lai - Baseband Receiver Design for Wireless MIMO-OFDM Communications, Wiley, 2012.

## License

This project is for educational purposes. Classical algorithm implementations are based on [1] with modifications for experimenting and organizing.

## 🤖 AI Assistance

Parts of this project were developed with the assistance of AI tools (e.g., ChatGPT,Claude) for:
- Code suggestions and debugging
- Concept explanations
- Documentation drafting

All outputs were reviewed and adapted as needed.
