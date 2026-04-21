function x_GI = build_frame(Xmod, Nfft, Ng, Nused, Nvc, Nframe)
% Concatenate Nframe OFDM symbols (mapping + IFFT + CP) into one vector.
    x_GI = zeros(1, Nframe*(Nfft+Ng));
    for k = 1:Nframe
        idx    = (k-1)*Nused + (1:Nused);
        Xk     = Xmod(idx);
        X_shift = map_to_subcarriers(Xk, Nfft, Nused, Nvc);
        x_sym  = ofdm_mod(X_shift, Ng);
        x_GI((k-1)*(Nfft+Ng) + (1:Nfft+Ng)) = x_sym;
    end
end