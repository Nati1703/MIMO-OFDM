function Y_freq = ofdm_demod(y_time, N, Ncp)
% OFDM demodulation: CP removal + FFT
    Y_freq = fft(y_time(Ncp+1 : Ncp+N));
end
