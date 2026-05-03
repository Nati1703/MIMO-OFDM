%% Channel estimation: static SISO OFDM
clear; close all; clc;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
addpath('../Channel_models');
addpath('../OFDM');
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);

addpath(repo_root);
if exist(fullfile(repo_root, 'channel'), 'dir')
    addpath(fullfile(repo_root, 'channel'));
end
if exist(fullfile(repo_root, 'Channel_models'), 'dir')
    addpath(fullfile(repo_root, 'Channel_models'));
end

% -------------------------------------------------------------------------
% Shared OFDM / modulation parameters
% -------------------------------------------------------------------------
Nfft  = 64;
Ng    = Nfft/4;
Nvc   = Nfft/4;
Nused = Nfft - Nvc;

Nbps = 2;
M    = 2^Nbps;

pilot_val = exp(1j*pi/4);              % constant unit-power pilot

% Channel parameters
Ts    = 1e-6;                           % sample period
fs    = 1 / Ts;
tau_d = 2 * Ts;
[pdp, lmax] = exp_pdp(tau_d, Ts);
L = lmax + 1;

fprintf('Nfft=%d, Ng=%d, Nused=%d, Nvc=%d, L=%d\n', Nfft, Ng, Nused, Nvc, L);

% -------------------------------------------------------------------------
% Sanity check: verify that the used-subcarrier indices recovered by
% probing map_to_subcarriers really are natural FFT bins, by checking
% that A * h reconstructs demap_from_subcarriers(fft(h, Nfft), ...) for
% an arbitrary length-L channel. If map/demap ever drift to inconsistent
% conventions (one shifted, one not), this fires.
% -------------------------------------------------------------------------
h_chk = randn(1, L) + 1j*randn(1, L);
Hfull_chk = fft([h_chk zeros(1, Nfft-L)], Nfft);
Hused_chk = demap_from_subcarriers(Hfull_chk, Nfft, Nused, Nvc);

probe_chk = 1:Nused;
probe_full_chk = map_to_subcarriers(probe_chk, Nfft, Nused, Nvc);
probe_int_chk = round(real(probe_full_chk));
k_used_chk = zeros(1, Nused);
for m = 1:Nused
    k_used_chk(m) = find(probe_int_chk == m, 1) - 1;
end
A_chk = exp(-1j * 2*pi * (k_used_chk(:) * (0:L-1)) / Nfft);
idx_err = max(abs(Hused_chk(:) - A_chk * h_chk(:)));
assert(idx_err < 1e-9, ...
    'Subcarrier indexing mismatch (err = %g). Check map_to_subcarriers / demap_from_subcarriers conventions.', idx_err);
fprintf('Subcarrier-indexing check passed (max err = %.2e)\n', idx_err);

%% Section 1: Single-shot static channel estimate
EbN0_dB_1 = 20;
sigma2_t = 1 / (Nfft * Nbps * 10^(EbN0_dB_1/10));   % time-domain noise variance
noise_var_fd = Nfft * sigma2_t;                     % FFT-bin noise variance

h = sqrt(pdp(:).') .* Rayleigh_model(L);
Hfull_true = fft([h zeros(1, Nfft-L)], Nfft);
Htrue_used = demap_from_subcarriers(Hfull_true, Nfft, Nused, Nvc);

% all-known pilot OFDM symbol
Xpilot_used = pilot_val * ones(1, Nused);
x_pilot = tx_ofdm_symbol(Xpilot_used, Nfft, Ng, Nused, Nvc);
y_pilot = apply_channel_symbol(x_pilot, h, sigma2_t);
Ypilot_used = rx_ofdm_symbol(y_pilot, Nfft, Ng, Nused, Nvc);

Hls    = ls_ce(Ypilot_used, Xpilot_used);
Hmmse  = mmse_ce(Hls, noise_var_fd, pdp, Nfft, Nvc);
Hdftls = dft_ce(Hls, L, Nfft, Nvc);

figure('Name', 'Static CE: single-shot', 'Position', [100 100 1200 450]);

subplot(1,2,1);
plot(abs(Htrue_used), 'k-', 'LineWidth', 1.8); hold on;
plot(abs(Hls), '--', 'LineWidth', 1.2);
plot(abs(Hmmse), '-.', 'LineWidth', 1.2);
plot(abs(Hdftls), ':', 'LineWidth', 1.8);
grid on;
xlabel('Used-subcarrier index');
ylabel('|H|');
title(sprintf('Static channel estimate @ %d dB', EbN0_dB_1));
legend('True', 'LS', 'MMSE', 'DFT-LS', 'Location', 'best');

subplot(1,2,2);
% Compute the L-tap LS estimate of the channel taps (this is the
% intermediate quantity that dft_ce solves for internally).
A_S1 = exp(-1j * 2*pi * (k_used_chk(:) * (0:L-1)) / Nfft);
h_est_ls = A_S1 \ Hls(:);
stem(0:L-1, abs(h),         'k',  'LineWidth', 1.8); hold on;
stem(0:L-1, abs(h_est_ls),  'd--', 'LineWidth', 1.5);
xlim([-0.5, L-0.5]);
grid on;
xlabel('Tap index');
ylabel('|h|');
title('Channel impulse response');
legend('True', 'L-tap LS estimate', 'Location', 'best');

%% Section 2: Static channel MSE vs SNR
EbN0_dB = 0:2:24;
Niter_mse = 400;

mse_ls    = zeros(size(EbN0_dB));
mse_mmse  = zeros(size(EbN0_dB));
mse_dftls = zeros(size(EbN0_dB));

for i = 1:numel(EbN0_dB)
    sigma2_t = 1 / (Nfft * Nbps * 10^(EbN0_dB(i)/10));
    noise_var_fd = Nfft * sigma2_t;

    acc_ls = 0;
    acc_mmse = 0;
    acc_dft = 0;

    for it = 1:Niter_mse
        h = sqrt(pdp(:).') .* Rayleigh_model(L);
        Htrue_used = demap_from_subcarriers(fft([h zeros(1, Nfft-L)], Nfft), Nfft, Nused, Nvc);

        x_pilot = tx_ofdm_symbol(Xpilot_used, Nfft, Ng, Nused, Nvc);
        y_pilot = apply_channel_symbol(x_pilot, h, sigma2_t);
        Ypilot_used = rx_ofdm_symbol(y_pilot, Nfft, Ng, Nused, Nvc);

        Hls    = ls_ce(Ypilot_used, Xpilot_used);
        Hmmse  = mmse_ce(Hls, noise_var_fd, pdp, Nfft, Nvc);
        Hdftls = dft_ce(Hls, L, Nfft, Nvc);

        acc_ls   = acc_ls   + mean(abs(Hls    - Htrue_used).^2);
        acc_mmse = acc_mmse + mean(abs(Hmmse  - Htrue_used).^2);
        acc_dft  = acc_dft  + mean(abs(Hdftls - Htrue_used).^2);
    end

    mse_ls(i)    = acc_ls   / Niter_mse;
    mse_mmse(i)  = acc_mmse / Niter_mse;
    mse_dftls(i) = acc_dft  / Niter_mse;

    fprintf('Static MSE done: Eb/N0 = %2d dB\n', EbN0_dB(i));
end

figure('Name', 'Static CE MSE', 'Position', [100 100 650 500]);
semilogy(EbN0_dB, mse_ls, 'o-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, mse_mmse, 's-', 'LineWidth', 1.5);
semilogy(EbN0_dB, mse_dftls, 'd-', 'LineWidth', 1.5);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Channel MSE');
title('Static SISO OFDM channel estimation');
legend('LS', 'MMSE', 'DFT-LS', 'Location', 'southwest');

%% Section 3: Static BER vs SNR (preamble estimate + one-tap EQ)
NdataSym_static = 6;
Target_neb = 600;
Niter_ber = 4000;

ber_ls    = zeros(size(EbN0_dB));
ber_mmse  = zeros(size(EbN0_dB));
ber_dftls = zeros(size(EbN0_dB));

for i = 1:numel(EbN0_dB)
    sigma2_t = 1 / (Nfft * Nbps * 10^(EbN0_dB(i)/10));
    noise_var_fd = Nfft * sigma2_t;

    Neb_ls = 0; Neb_mmse = 0; Neb_dft = 0;
    Ntb = 0;

    for it = 1:Niter_ber
        h = sqrt(pdp(:).') .* Rayleigh_model(L);

        % Preamble
        x_pilot = tx_ofdm_symbol(Xpilot_used, Nfft, Ng, Nused, Nvc);
        y_pilot = apply_channel_symbol(x_pilot, h, sigma2_t);
        Ypilot_used = rx_ofdm_symbol(y_pilot, Nfft, Ng, Nused, Nvc);

        Hls    = ls_ce(Ypilot_used, Xpilot_used);
        Hmmse  = mmse_ce(Hls, noise_var_fd, pdp, Nfft, Nvc);
        Hdftls = dft_ce(Hls, L, Nfft, Nvc);

        for s = 1:NdataSym_static
            bits_tx = randi([0 1], 1, Nused * Nbps);
            Xused   = qam_mod(bits_tx, M);

            x_sym = tx_ofdm_symbol(Xused, Nfft, Ng, Nused, Nvc);
            y_sym = apply_channel_symbol(x_sym, h, sigma2_t);
            Yused = rx_ofdm_symbol(y_sym, Nfft, Ng, Nused, Nvc);

            Xhat_ls   = one_tap_eq(Yused, Hls);
            Xhat_mmse = one_tap_eq(Yused, Hmmse); %% all are ZF equalization approaches
            Xhat_dft  = one_tap_eq(Yused, Hdftls);

            bits_ls   = qam_demod(Xhat_ls, M);
            bits_mmse = qam_demod(Xhat_mmse, M);
            bits_dft  = qam_demod(Xhat_dft, M);

            Neb_ls   = Neb_ls   + sum(bits_tx ~= bits_ls);
            Neb_mmse = Neb_mmse + sum(bits_tx ~= bits_mmse);
            Neb_dft  = Neb_dft  + sum(bits_tx ~= bits_dft);
            Ntb = Ntb + numel(bits_tx);
        end

        if max([Neb_ls Neb_mmse Neb_dft]) > Target_neb
            break;
        end
    end

    ber_ls(i)    = Neb_ls   / Ntb;
    ber_mmse(i)  = Neb_mmse / Ntb;
    ber_dftls(i) = Neb_dft  / Ntb;

    fprintf('Static BER done: Eb/N0 = %2d dB\n', EbN0_dB(i));
end

figure('Name', 'Static BER', 'Position', [100 100 650 500]);
semilogy(EbN0_dB, ber_ls, 'o-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, ber_mmse, 's-', 'LineWidth', 1.5);
semilogy(EbN0_dB, ber_dftls, 'd-', 'LineWidth', 1.5);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('BER');
title('Static SISO OFDM BER');
legend('LS', 'MMSE', 'DFT-LS', 'Location', 'southwest');
ylim([1e-5 1]);

%% Local functions
function x_sym = tx_ofdm_symbol(Xused, Nfft, Ng, Nused, Nvc)
    Xshift = map_to_subcarriers(Xused, Nfft, Nused, Nvc);
    x_sym = ofdm_mod(Xshift, Ng);
end

function Yused = rx_ofdm_symbol(y_sym, Nfft, Ng, Nused, Nvc)
    Yshift = ofdm_demod(y_sym, Nfft, Ng);
    Yused = demap_from_subcarriers(Yshift, Nfft, Nused, Nvc);
end

function y = apply_channel_symbol(x, h, sigma2_t)
    y = conv(x, h);
    y = y(1:numel(x));
    if sigma2_t > 0
        y = y + sqrt(sigma2_t/2) * (randn(size(y)) + 1j*randn(size(y)));
    end
end


%{
The same channel-estimation ideas that we saw here naturally extend to
MIMO-OFDM. Instead of a scalar per subcarrier they become channel matrices.
Revisit Weekly guidance from Ivlrac's notes (g04) for some aspects of why
orthogonal pilot sequences are necessary. All three(LS, MMSE, and DFT-LS)
extend naturally.
%}