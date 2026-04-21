function symbols = qam_mod(bits, M)
% QAM modulation: bits → unit-power complex symbols (gray coded)
    Nbps = log2(M);
    N_sym = length(bits) / Nbps;
    bit_groups = reshape(bits, Nbps, N_sym).';
    int_vals = bi2de(bit_groups, 'left-msb'); % matlab by default right-msb
    symbols = qammod(int_vals, M, 'gray', 'UnitAveragePower', true).';
end