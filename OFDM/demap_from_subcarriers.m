function Xmod = demap_from_subcarriers(X_shift, Nfft, Nused, Nvc)
% Inverse of map_to_subcarriers.
    Xmod = [X_shift(Nfft-Nused/2+1:Nfft), X_shift(2:Nused/2+1)];
end
