function H = kronecker_MIMO(Nr, Nt, R_tx, R_rx, varargin)
% Generate spatially correlated MIMO channel (Kronecker model)
% H = R_rx^(1/2) * H_iid * R_tx^(1/2)
%
% default: Rayleigh; pass in K_dB for Rician
%
% Inputs:
%   Nr, Nt  - number of Rx/Tx antennas
%   R_tx    - Nt x Nt transmit correlation matrix
%   R_rx    - Nr x Nr receive correlation matrix
%   K_dB    - (optional) Rician K-factor in dB. If omitted, Rayleigh is used.

    if isempty(varargin)
        H_iid = Rayleigh_MIMO(Nr, Nt);
    else
        K_dB = varargin{1};
        H_iid = Rician_MIMO(Nr, Nt, K_dB);
    end
    
    H = sqrtm(R_rx) * H_iid * sqrtm(R_tx);
end