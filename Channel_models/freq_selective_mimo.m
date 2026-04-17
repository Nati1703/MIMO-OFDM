function [H_freq, h_taps] = freq_selective_mimo(Nr, Nt, N, tau_d, Ts, R_tx, R_rx)
% Generate frequency-selective MIMO channel
%
% Inputs:
%   Nr, Nt  - number of Rx/Tx antennas
%   N       - FFT size (number of subcarriers)
%   tau_d   - RMS delay spread [seconds]
%   Ts      - Sampling period [seconds]
%   R_tx    - Tx correlation matrix (eye(Nt) for uncorrelated)
%   R_rx    - Rx correlation matrix (eye(Nr) for uncorrelated)
%
% Outputs:
%   H_freq  - Nr x Nt x N channel in frequency domain
%   h_taps  - Nr x Nt x L channel taps in delay domain
%
% Steps:
%   1. Generate exponential PDP (determines L automatically)
%   2. At each tap: generate correlated MIMO matrix, scale by sqrt(PDP)
%   3. FFT across taps → H[k] per subcarrier
 
    [pdp, lmax] = exp_pdp(tau_d, Ts);
    L = lmax + 1;
    
    h_taps = zeros(Nr, Nt, L);
    for l = 1:L
        h_taps(:,:,l) = sqrt(pdp(l)) * kronecker_MIMO(Nr, Nt, R_tx, R_rx);
    end
    
    H_freq = zeros(Nr, Nt, N);
    for rx = 1:Nr
        for tx = 1:Nt
            h_vec = squeeze(h_taps(rx, tx, :));
            H_freq(rx, tx, :) = fft(h_vec, N);
        end
    end
end
 