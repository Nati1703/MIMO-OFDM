function y = doppler_spectrum(fd, Nfft)
% Clarke/Jakes Doppler power spectrum (Cho Program 2.13)
% Input:  fd   = maximum Doppler frequency [Hz]
%         Nfft = number of frequency domain points
% Output: y    = Doppler spectrum (1 x Nfft)
    df = 2*fd / Nfft;
    f(1) = 0;
    y(1) = 1.5 / (pi*fd);
    for i = 2:Nfft/2
        f(i) = (i-1) * df;
        y([i Nfft-i+2]) = 1.5 / (pi*fd*sqrt(1-(f(i)/fd)^2));
    end
    nFitPoints = 3;
    kk = (Nfft/2 - nFitPoints) : Nfft/2;
    polyFreq = polyfit(f(kk), y(kk), nFitPoints);
    y(Nfft/2 + 1) = polyval(polyFreq, f(Nfft/2) + df);
end