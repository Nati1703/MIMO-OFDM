function bits = qam_demod(symbols, M)
% QAM demodulation: complex symbols → bits (hard decision, gray coded)
    Nbps = log2(M);
    int_vals = qamdemod(symbols.', M, 'gray', 'UnitAveragePower', true);
    bit_groups = de2bi(int_vals, Nbps, 'left-msb');
    bits = reshape(bit_groups.', 1, []);
end