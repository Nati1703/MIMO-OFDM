function [h, Nfft, Nifft, doppler_coeff] = fwgn_model(fm, fs, N)
% FWGN (Clarke/Gans) model (Cho Program 2.12)
% Generates time-correlated Rayleigh fading
%
% Inputs:
%   fm  = maximum Doppler frequency [Hz]
%   fs  = sampling frequency [Hz]
%   N   = number of output samples
%
% Output:
%   h   = complex fading sequence (1 x N), unit average power
    Nfft = 2^max(3, nextpow2(2*fm/fs*N));
    Nifft = ceil(Nfft * fs / (2*fm));
    GI = randn(1, Nfft);
    GQ = randn(1, Nfft);
    CGI = fft(GI);
    CGQ = fft(GQ);
    doppler_coeff = doppler_spectrum(fm, Nfft);
    f_CGI = CGI .* sqrt(doppler_coeff);
    f_CGQ = CGQ .* sqrt(doppler_coeff);
    Filtered_CGI = [f_CGI(1:Nfft/2) zeros(1, Nifft-Nfft) f_CGI(Nfft/2+1:Nfft)];
    Filtered_CGQ = [f_CGQ(1:Nfft/2) zeros(1, Nifft-Nfft) f_CGQ(Nfft/2+1:Nfft)];
    hI = ifft(Filtered_CGI);
    hQ = ifft(Filtered_CGQ);
    rayEnvelope = sqrt(abs(hI).^2 + abs(hQ).^2);
    rayRMS = sqrt(mean(rayEnvelope(1:N) .* rayEnvelope(1:N)));
    h = complex(real(hI(1:N)), -real(hQ(1:N))) / rayRMS;
end