%% MIMO Detection: Verification
clear; close all; clc;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
addpath('../Channel_models');
addpath('../OFDM');

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end
repo_root = fileparts(this_dir);

addpath(this_dir);     
addpath(repo_root);

if exist(fullfile(repo_root, 'channel'), 'dir')
    addpath(fullfile(repo_root, 'channel'));
end
if exist(fullfile(repo_root, 'Channel_models'), 'dir')
    addpath(fullfile(repo_root, 'Channel_models'));
end
if exist(fullfile(repo_root, 'OFDM'), 'dir')
    addpath(fullfile(repo_root, 'OFDM'));
end
%% ========================================================================
%  SECTION 1: BER curves — ZF vs MMSE vs OSIC
%  ========================================================================

Nr = 2;
Nt = 2;
Nbps = 2;
M = 2^Nbps;

all_bits = de2bi(0:M-1, Nbps, 'left-msb');
constellation = qam_mod(reshape(all_bits.', 1, []), M).';

SNR_dB = 0:2:30;
N_frames = 500000;

BER_zf = zeros(size(SNR_dB));
BER_mmse = zeros(size(SNR_dB));
BER_osic = zeros(size(SNR_dB));

fprintf('=== SECTION 2: BER simulation ===\n');
for idx_snr = 1:length(SNR_dB)
    snr_lin = 10^(SNR_dB(idx_snr)/10);
    noise_var = Nt / snr_lin;  % per-receive-antenna noise variance

    err_zf = 0;
    err_mmse = 0;
    err_osic = 0;
    n_bits = 0;

    for frame = 1:N_frames
        % Channel model imported from the project utilities
        H = Rayleigh_MIMO(Nr, Nt);

        % Transmit random QPSK from each antenna using the project modem
        bits_tx = randi([0 1], 1, Nt * Nbps);
        x = qam_mod(bits_tx, M).';

        % Receive
        n = sqrt(noise_var/2) * (randn(Nr, 1) + 1j * randn(Nr, 1));
        y = H * x + n;

        % ZF
        x_zf = zf_detector(y, H);
        bits_zf = qam_demod(x_zf.', M);
        err_zf = err_zf + sum(bits_tx ~= bits_zf);

        % MMSE
        x_mmse = mmse_detector(y, H, noise_var);
        bits_mmse = qam_demod(x_mmse.', M);
        err_mmse = err_mmse + sum(bits_tx ~= bits_mmse);

        % OSIC already returns hard constellation decisions
        x_osic = osic_detector(y, H, noise_var, constellation);
        bits_osic = qam_demod(x_osic.', M);
        err_osic = err_osic + sum(bits_tx ~= bits_osic);

        n_bits = n_bits + numel(bits_tx);
    end

    BER_zf(idx_snr) = err_zf / n_bits;
    BER_mmse(idx_snr) = err_mmse / n_bits;
    BER_osic(idx_snr) = err_osic / n_bits;

    fprintf('SNR=%2d dB: ZF=%.3e  MMSE=%.3e  OSIC=%.3e\n', ...
        SNR_dB(idx_snr), BER_zf(idx_snr), BER_mmse(idx_snr), BER_osic(idx_snr));
end
%%
figure('Name', 'Detection Comparison', 'Position', [100 100 700 500]);
semilogy(SNR_dB, BER_zf, 'b-o', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, BER_mmse, 'r-s', 'LineWidth', 1.5);
semilogy(SNR_dB, BER_osic, 'g-^', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('BER');
title('2x2 MIMO Detection: ZF vs MMSE vs OSIC');
legend('ZF', 'MMSE', 'OSIC', 'Location', 'southwest');
grid on;
ylim([1e-5 1]);

