# MIMO-OFDM Link-Level Simulator

A modular link-level simulator for MIMO-OFDM systems, built as a living
companion to graduate coursework at the Technical University of Munich.
The goal is a plug-and-play platform: each module can be swapped, extended,
or replaced as new algorithms are studied in lectures.

Classical algorithms follow Cho et al. [1]. The project grows alongside
ongoing courses — see **Planned Extensions** below.

---

## System Architecture

```
TRANSMITTER
─────────────────────────────────────────────────────────
  Bits → QAM Mapper → S/P
                        ↓
      Antenna 1: IFFT → Add CP → Transmit
      Antenna 2: IFFT → Add CP → Transmit

CHANNEL
─────────────────────────────────────────────────────────
  Frequency-selective MIMO channel (Kronecker correlation)
  y[k] = H[k] · x[k] + n[k]   per subcarrier k

RECEIVER
─────────────────────────────────────────────────────────
  Each Rx antenna: Remove CP → FFT → Y_rx[k]
                                        ↓
         ┌──────────────────────────────────────────┐
         │  Per subcarrier k = 0, …, N-1:           │
         │  1. Channel Estimation  (LS / MMSE / DFT)│
         │  2. MIMO Detection      (ZF / MMSE / OSIC)│
         └──────────────────────────────────────────┘
                        ↓
               QAM Demod → Bits

Assumptions: *perfect synchronisation*, uncoded BER.
Synch issues(STO,CFO,SCO) are part of broader review but not in simulations.
```

---

## Repository Structure

```
mimo-ofdm/
│
├── README.md
│
├── Assets/
│   ├── OFDM_AWGN_16QAM.png       % BER curve: OFDM vs theory, 16-QAM
│   ├── OFDM_AWGN_256QAM.png      % BER curve: OFDM vs theory, 256-QAM
│   ├── ZP_CP.png                  % CP vs ZP MSE comparison
│   └── ZP_vs_CP.png               % CP vs ZP spectrum leakage
│
├── Preliminary_DFT_Review/
│   └── DFT_review.mlx             % DFT fundamentals: bins, orthogonality,
│                                  %   circular convolution, CP, windowing
│
├── Channel_models/
│   ├── Channels_verification.mlx  % Full statistical verification (all models)
│   ├── Rayleigh_model.m           % SISO Rayleigh fading samples
│   ├── Rayleigh_MIMO.m            % i.i.d. Rayleigh MIMO matrix
│   ├── Rician_model.m             % SISO Rician fading (K-factor in dB)
│   ├── Rician_MIMO.m              % Rician MIMO (configurable K)
│   ├── kronecker_MIMO.m           % Spatially correlated MIMO (Kronecker)
│   ├── exp_pdp.m                  % Exponential PDP (physical units: τ_d, T_s)
│   ├── freq_selective_mimo.m      % Static freq-selective MIMO channel
│   ├── doppler_spectrum.m         % Clarke/Jakes Doppler spectrum
│   ├── fwgn_model.m               % Time-correlated fading per tap (FWGN)
│   └── time_varying_mimo.m        % Doubly-dispersive MIMO-OFDM channel
│
├── OFDM/
│   ├── OFDM_review.m              % Verification: QAM, loopback, BER, CP vs ZP
│   ├── qam_mod.m                  % QAM modulation (Gray-coded, unit power)
│   ├── qam_demod.m                % QAM demodulation (hard decision)
│   ├── ofdm_mod.m                 % IFFT + cyclic prefix insertion
│   ├── ofdm_demod.m               % CP removal + FFT
│   ├── map_to_subcarriers.m       % Data symbols → FFT bins (guard band)
│   ├── demap_from_subcarriers.m   % FFT bins → data symbols
│   ├── build_frame.m              % Build multi-symbol OFDM frame
│   ├── extract_frame.m            % Extract symbols from received frame
│   └── ber_QAM.m                  % Analytical BER for M-QAM (AWGN/Rayleigh)
│
├── Channel_estimation/
│   ├── channel_estimation.m       % Verification: MSE and BER for all estimators
│   ├── ls_ce.m                    % LS estimation: Ĥ = Y / X
│   ├── mmse_ce.m                  % LMMSE estimation (R_HH from PDP)
│   ├── dft_ce.m                   % DFT-based denoising: IFFT → keep L → FFT
│   └── one_tap_eq.m               % Per-subcarrier ZF equalisation: Y / Ĥ
│
└── Channel_detection/
    ├── detection_verification.m   % Verification: BER for ZF, MMSE, OSIC
    ├── zf_detector.m              % Zero-forcing: (H'H)⁻¹H'y
    ├── mmse_detector.m            % MMSE: (H'H + σ²I)⁻¹H'y
    └── osic_detector.m            % OSIC with SINR-based ordering
```

---

## What Each Verification Script Produces

### `Preliminary_DFT_Review/DFT_review.mlx`
Builds intuition for the DSP concepts underlying OFDM:
- DFT/IDFT roundtrip, frequency axis construction, bin spacing
- Orthogonality of complex exponentials (integer vs non-integer frequency)
- Linear vs circular convolution
- Cyclic prefix: Y[k] = H[k]X[k] with CP, fails without
- Spectral leakage and windowing (rectangular vs raised cosine)
- Channel frequency response from tap coefficients

### `Channel_models/Channels_verification.mlx`
Statistical tests for every channel model:
- **Section 1:** Rayleigh envelope PDF, I/Q scatter, phase uniformity, power distribution
- **Section 1b:** Rician PDF for K = 0, 1, 5, 20 — I/Q clouds shift toward LOS
- **Section 2:** Exponential PDP tap powers vs theory, physical delay axis (ns)
- **Section 3:** Kronecker — singular value distributions vs ρ, correlation matrix error
- **Section 4:** Frequency-selective — impulse response, frequency response, Parseval check
- **Section 5:** Monte Carlo — average subcarrier power ≈ 1, frequency correlation r_f[Δk]
- **Section 5b:** FWGN — Rayleigh envelope, J₀ autocorrelation, coherence time table
- **Section 5c:** Doubly-dispersive — time-frequency waterfall plot H[k,t]

### `OFDM/OFDM_review.m`
- **Section 1:** QAM for M = 4, 16, 64 — constellations with Gray-coded bit labels
- **Section 2:** OFDM mod/demod loopback — roundtrip error ≈ 0, CP copy verification
- **Section 3:** AWGN BER — 16-QAM simulation matches analytic curve exactly
- **Section 4:** CP vs zero-padding — per-subcarrier MSE and single-tone ICI leakage

### `Channel_estimation/channel_estimation.m`
Static SISO-OFDM with full preamble (block-type pilots):
- **Section 1:** Single-shot — true vs LS vs MMSE vs DFT-LS channel magnitude + tap estimates
- **Section 2:** MSE vs Eb/N₀ — MMSE and DFT-LS both beat LS at low SNR
- **Section 3:** BER vs Eb/N₀ — end-to-end with one-tap ZF equalisation

### `Channel_detection/detection_verification.m`
Flat Rayleigh 2×2 MIMO, QPSK:
- **Section 1:** BER curves — ZF < MMSE < OSIC ordering confirmed

---



## Planned Extensions

The simulator is designed to grow alongside my current courses:

### Modern Channel Codes *(TUM)*
Adds a soft-output receiver and channel coding layer:

```
Current:   Estimation → Detection → Bits            (uncoded, one-shot)

Planned:   Estimation ⇄ Detection ⇄ Decoder         (iterative / turbo)
```

Specific planned additions (incremental — standard link-sim building blocks):
- **Soft-output detection:** replace hard QAM decisions with **per-bit LLRs** (log-likelihood ratios) after MIMO equalisation — same interface a channel decoder expects as input.
- **Channel decoding:** feed those LLRs into an **LDPC** or **Turbo** decoder (MAP/BCJR- or belief-propagation–style APP decoding, depending on the code), and compare coded vs uncoded BER.
- **Optional outer loop:** a few **detector ⇄ decoder** iterations (extrinsic information exchange) — scoped to what stays tractable in this simulator.
---

## Dependencies

- **MATLAB** R2020a or later (Communications Toolbox for `qammod`/`qamdemod`)

---

## References

1. Y.S. Cho et al., *MIMO-OFDM Wireless Communications with MATLAB*, Wiley, 2010.
2. T.-D. Chiueh, P.-Y. Tsai, *Baseband Receiver Design for Wireless MIMO-OFDM*, Wiley, 2012.

---

## AI Assistance

Parts of this project were developed with the assistance of Chatgpt
for code suggestions, concept explanations, and documentation drafting.
All outputs were reviewed, tested, and adapted.

---
