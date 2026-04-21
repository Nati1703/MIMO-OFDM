function ber = ber_QAM_clean(EbN0dB, M_qam, channel_type)
% Analytical BER for square M-QAM (Gray coded)
% Expressions based on (Cho Eq. 4.25 / 4.26)

    L = sqrt(M_qam);                 % levels per I/Q axis
    k = log2(M_qam);                 % bits per symbol
    gamma_b = 10.^(EbN0dB/10);       % Eb/N0 (linear)

    if lower(channel_type(1)) == 'a'
        % ===== AWGN channel =====
        a = 4/k * (1 - 1/L);
        x = 3*k/(M_qam - 1) .* gamma_b;
        % Q(x) = 0.5 * erfc(x / sqrt(2))
        ber = a * 0.5 * erfc( sqrt(x/2) );
    else
        % ===== Rayleigh fading =====
        % Average BER over Rayleigh fading
        x = 3*k/(M_qam - 1) .* gamma_b;
        ber = (2/k)*(1 - 1/L) .* ...
              (1 - sqrt( x ./ (x + 2) ));
    end
end