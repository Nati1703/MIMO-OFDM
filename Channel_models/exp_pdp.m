function [pdp, lmax] = exp_pdp(tau_d, Ts, A_dB, norm_flag)
% Generate normalized exponential power delay profile 
%
% Inputs:
%   tau_d     : RMS delay spread [seconds]
%   Ts        : Sampling period [seconds]
%   A_dB      : Smallest noticeable power in dB (default: -20)
%   norm_flag : Normalize total power to 1 (default: true)
%
% Outputs:
%   pdp       : Power delay profile vector
%   lmax      : Number of taps - 1

 
    if nargin < 4, norm_flag = 1; end
    if nargin < 3, A_dB = -20; end
    
    sigma_tau = tau_d;
    A = 10^(A_dB/10);
    lmax = ceil(-tau_d * log(A) / Ts);  % number of taps
    
    l = 0:lmax;
    
    if norm_flag
        % Normalize so total discrete power sums to 1
        p0 = (1 - exp(-Ts/sigma_tau)) / (1 - exp(-(lmax+1)*Ts/sigma_tau));
    else
        p0 = 1/sigma_tau;
    end
    
    pdp = p0 * exp(-l * Ts / sigma_tau);
end