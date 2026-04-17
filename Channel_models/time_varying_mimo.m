function [H_freq_tv, h_taps_tv] = time_varying_mimo(Nr, Nt, N, T, tau_d, Ts, fm, fs_ofdm, R_tx, R_rx)
% Generate time-varying frequency-selective MIMO channel
%
% Inputs:
%   Nr, Nt    - number of Rx/Tx antennas
%   N         - FFT size
%   T         - number of OFDM symbols (time snapshots)
%   tau_d     - RMS delay spread [seconds]
%   Ts        - sampling period [seconds]
%   fm        - maximum Doppler frequency [Hz]
%   fs_ofdm   - OFDM symbol rate [Hz]
%   R_tx, R_rx - correlation matrices
%
% Outputs:
%   H_freq_tv - Nr x Nt x N x T channel
%   h_taps_tv - Nr x Nt x L x T channel taps
    [pdp, lmax] = exp_pdp(tau_d, Ts);
    L = lmax + 1;
    h_taps_tv = zeros(Nr, Nt, L, T);
    for rx = 1:Nr
        for tx = 1:Nt
            for l = 1:L
                h_time = fwgn_model(fm, fs_ofdm, T);
                h_taps_tv(rx, tx, l, :) = sqrt(pdp(l)) * h_time;
            end
        end
    end
    H_freq_tv = zeros(Nr, Nt, N, T);
    for t = 1:T
        for rx = 1:Nr
            for tx = 1:Nt
                h_vec = squeeze(h_taps_tv(rx, tx, :, t));
                H_freq_tv(rx, tx, :, t) = fft(h_vec, N);
            end
        end
    end
end