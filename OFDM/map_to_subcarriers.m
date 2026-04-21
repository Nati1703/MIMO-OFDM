function X_shift = map_to_subcarriers(Xmod, Nfft, Nused, Nvc)
% Cho-style mapping: DC null, positive bins, VC nulls, negative bins.
    X_shift = zeros(1, Nfft);
    X_shift(2:Nused/2+1) = Xmod(Nused/2+1:Nused);
    X_shift(Nfft-Nused/2+1:Nfft)    = Xmod(1:Nused/2);

end