function x_time = ofdm_mod(X_freq, Ncp) % assumes X_freq already frequency aligned correctly
% OFDM modulation: IFFT + cyclic prefix
    x = ifft(X_freq);
    x_time = [x(end-Ncp+1:end), x];
end