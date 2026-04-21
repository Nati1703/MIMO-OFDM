function Xmod_r = extract_frame(y_GI, Nfft, Ng, Nused, Nvc, Nframe)
% Inverse of build_frame for a noiseless/AWGN channel (no equalisation).
    Xmod_r = zeros(1, Nframe*Nused);
    for k = 1:Nframe
        y_sym  = y_GI((k-1)*(Nfft+Ng) + (1:Nfft+Ng));
        Y      = ofdm_demod(y_sym, Nfft, Ng);
        Xmod_r((k-1)*Nused + (1:Nused)) = demap_from_subcarriers(Y, Nfft, Nused, Nvc);
    end
end