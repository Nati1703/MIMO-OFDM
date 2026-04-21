%% OFDM review: QAM sanity, modem loopback, AWGN BER vs theory, CP vs ZP
clear; close all; clc;
addpath('../Channel_models');


%% Section 1: QAM sanity
N_test = 10000;

for M = [4, 16, 64]
    Nbps = log2(M);
    bits_in  = randi([0 1], 1, N_test * Nbps);
    symbols  = qam_mod(bits_in, M);
    avg_pow  = mean(abs(symbols).^2);
    bits_out = qam_demod(symbols, M);
    n_err    = sum(bits_in ~= bits_out);
    fprintf('%d-QAM: %d symbols, avg power = %.4f (target 1), bit errors = %d\n', ...
        M, length(symbols), avg_pow, n_err);
end

figure('Name', 'Section 1: QAM constellations', 'Position', [100 100 1200 350]);
M_list = [4, 16, 64];
for idx = 1:3
    M    = M_list(idx);
    Nbps = log2(M);
    all_bits    = de2bi(0:M-1, Nbps, 'left-msb');
    all_bits_v  = reshape(all_bits.', 1, []);
    all_symbols = qam_mod(all_bits_v, M);
    subplot(1,3,idx);
    plot(real(all_symbols), imag(all_symbols), 'b.', 'MarkerSize', 15);
    axis equal; grid on;
    title(sprintf('%d-QAM (%d bits/symbol)', M, Nbps));
    xlabel('In-phase'); ylabel('Quadrature');
    if M <= 16
        for s = 1:M
            lbl = strrep(num2str(all_bits(s,:)), ' ', '');
            text(real(all_symbols(s))+0.05, imag(all_symbols(s))+0.05, lbl, 'FontSize', 7);
        end
    end
end
fprintf('\n');

%% Section 2: OFDM mod/demod loopback
Nfft = 64; Ng = Nfft/4; Nvc = Nfft/4;
Nused = Nfft - Nvc;
M = 4; Nbps = log2(M);

bits    = randi([0 1], 1, Nused * Nbps);
Xmod    = qam_mod(bits, M);
X_shift = map_to_subcarriers(Xmod, Nfft, Nused, Nvc);

x_time  = ofdm_mod(X_shift, Ng);
Y_shift = ofdm_demod(x_time, Nfft, Ng);

Ymod           = demap_from_subcarriers(Y_shift, Nfft, Nused, Nvc);
bits_recovered = qam_demod(Ymod, M);

fprintf('=== Section 2: OFDM mod/demod ===\n');
fprintf('Nfft=%d  Ng=%d  Nused=%d  Nvc=%d\n', Nfft, Ng, Nused, Nvc);
fprintf('Roundtrip error = %.2e,  bit errors = %d\n', ...
    norm(Xmod - Ymod), sum(bits ~= bits_recovered));

figure('Name', 'Section 2: OFDM signal', 'Position', [100 100 1200 600]);
subplot(2,2,1);
stem(0:Nfft-1, abs(X_shift), 'filled', 'MarkerSize', 3);
xlabel('FFT bin k'); ylabel('|X[k]|'); title('Subcarrier mapping'); grid on;

subplot(2,2,2);
plot(0:Nfft+Ng-1, real(x_time), 'b-', 'LineWidth', 1);
xlabel('Sample n'); ylabel('Re\{x[n]\}');
title('Time-domain OFDM symbol with CP'); hold on;
xline(Ng, 'r--', 'CP end', 'LineWidth', 1.5); grid on;

subplot(2,2,3);
plot(0:Ng-1, real(x_time(1:Ng)), 'r-',  'LineWidth', 2); hold on;
plot(0:Ng-1, real(x_time(Nfft+1:Nfft+Ng)), 'b--', 'LineWidth', 1.5);
title('CP = copy of last Ng samples'); xlabel('Sample');
legend('First Ng (CP)', 'Last Ng (data)'); grid on;

subplot(2,2,4);
stem(abs(Xmod - Ymod), 'filled', 'MarkerSize', 3);
xlabel('Data subcarrier'); ylabel('Error');
title('Roundtrip error per subcarrier'); grid on;
fprintf('\n');

%% Section 3: AWGN BER vs theory
% Noise calibration (MATLAB ifft/fft, unit-power QAM on data bins):
%   E[|X[k]|^2]=1 on used subcarriers  =>  Es=1 at the FFT output.
%   Time samples w[n] ~ CN(0,sigma2) i.i.d.  =>  FFT bin noise variance Nfft*sigma2.
%   Target Es/N0 = Nbps*Eb/N0  =>  sigma2 = 1/(Nfft*Nbps*10^(EbN0/10)).


Nfft = 64; Ng = Nfft/4; Nvc = Nfft/4;
Nused = Nfft - Nvc;
Nbps = 4; M = 2^Nbps;      % 16-QAM
Nframe = 3; N_iter = 5000; Target_neb = 500;
EbN0_dB = 0:2:24;
BER_sim = zeros(size(EbN0_dB));

for i = 1:length(EbN0_dB)
    sigma2 = 1 / (Nfft * Nbps * 10^(EbN0_dB(i)/10));
    Neb = 0; Ntb = 0;
    for m = 1:N_iter
        % Tx
        bits_tx = randi([0 1], 1, Nused*Nframe*Nbps);
        Xmod    = qam_mod(bits_tx, M);
        x_GI    = build_frame(Xmod, Nfft, Ng, Nused, Nvc, Nframe);

        % AWGN
        noise = sqrt(sigma2/2) * (randn(size(x_GI)) + 1j*randn(size(x_GI)));
        y_GI  = x_GI + noise;

        % Rx
        Ymod = extract_frame(y_GI, Nfft, Ng, Nused, Nvc, Nframe);
        bits_rx = qam_demod(Ymod, M);

        Neb = Neb + sum(bits_tx ~= bits_rx);
        Ntb = Ntb + length(bits_tx);
        if Neb > Target_neb, break; end
    end
    BER_sim(i) = Neb / Ntb;
    fprintf('Eb/N0 = %2d dB:  BER = %.3e\n', EbN0_dB(i), BER_sim(i));
end

EbN0_fine = 0:0.5:30;
figure('Name', 'Section 3: AWGN BER', 'Position', [100 100 600 500]);
semilogy(EbN0_fine, ber_QAM(EbN0_fine, M, 'AWGN'), 'r-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, BER_sim, 'b-', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title(sprintf('OFDM in AWGN, %d-QAM', M));
legend('AWGN analytic', 'OFDM simulation', 'Location', 'southwest');
grid on; ylim([1e-6 1]);
fprintf('\n');

%% Section 4: CP vs ZP -- orthogonality and ICI under multipath
% Three receivers share the same Tx data.
%   CP       : standard, uses ofdm_mod / ofdm_demod.
%   ZP-naive : Tx appends Ng zeros; Rx FFTs the first Nfft samples directly.
%   ZP+OA    : same Tx; Rx folds the last Ng received samples back onto the
%              first Ng samples before the FFT.
%              Intuition: ZP gives linear convolution; overlap-add restores
%              the circular structure as long as channel length <= Ng+1.

Nfft = 64; Ng = Nfft/4; Nvc = Nfft/4; Nused = Nfft - Nvc;
M = 4; Nbps = log2(M);

% Channel: short exponential-PDP FIR
Ts = 1; 
tau_d = 2*Ts;
[pdp, lmax] = exp_pdp(tau_d, Ts);
L = lmax + 1;


assert(L <= Ng + 1, 'Channel length %d exceeds Ng+1=%d', L, Ng+1);

%% --- Experiment 4a: MSE per data subcarrier, averaged over realizations ---
N_real = 1000;
mse_cp  = zeros(1, Nused);
mse_zpn = zeros(1, Nused);
mse_zpo = zeros(1, Nused);

for r = 1:N_real
    % Channel realization
    h = zeros(1, L);
    for l = 1:L
        h(l) = sqrt(pdp(l)) * Rayleigh_MIMO(1,1);
    end

    H    = fft(h, Nfft);
    Hdat = H_data_bins(H, Nfft, Nused, Nvc);

    % Tx
    bits    = randi([0 1], 1, Nused*Nbps);
    Xmod    = qam_mod(bits, M);
    X_shift = map_to_subcarriers(Xmod, Nfft, Nused, Nvc);

    % CP symbol
    x_cp = ofdm_mod(X_shift, Ng);

    % ZP symbol: append zeros at the end
    x_core = ifft(X_shift);
    x_zp   = [x_core, zeros(1, Ng)];

    % Channel
    y_cp = conv(x_cp, h);
    y_cp = y_cp(1:Nfft+Ng);

    y_zp = conv(x_zp, h);
    y_zp = y_zp(1:Nfft+Ng);

    % Rx
    % CP: drop prefix, FFT
    Y_cp = fft(y_cp(Ng+1 : Ng+Nfft));

    % ZP-naive: take first Nfft samples directly
    Y_zpn = fft(y_zp(1:Nfft));

    % ZP + overlap-add
    r_oa = y_zp(1:Nfft);
    r_oa(1:Ng) = r_oa(1:Ng) + y_zp(Nfft+1 : Nfft+Ng);
    Y_zpo = fft(r_oa);

    % Per-data-bin error after one-tap equalization
    Xcp  = demap_from_subcarriers(Y_cp,  Nfft, Nused, Nvc) ./ Hdat;
    Xzpn = demap_from_subcarriers(Y_zpn, Nfft, Nused, Nvc) ./ Hdat;
    Xzpo = demap_from_subcarriers(Y_zpo, Nfft, Nused, Nvc) ./ Hdat;

    mse_cp  = mse_cp  + abs(Xcp  - Xmod).^2;
    mse_zpn = mse_zpn + abs(Xzpn - Xmod).^2;
    mse_zpo = mse_zpo + abs(Xzpo - Xmod).^2;
end

mse_cp  = mse_cp  / N_real;
mse_zpn = mse_zpn / N_real;
mse_zpo = mse_zpo / N_real;

fprintf('=== Section 4: CP vs ZP (noiseless multipath) ===\n');
fprintf('mean MSE over data subcarriers\n');
fprintf('  CP       : %.3e\n', mean(mse_cp));
fprintf('  ZP-naive : %.3e\n', mean(mse_zpn));
fprintf('  ZP+OA    : %.3e\n', mean(mse_zpo));

%% --- Experiment 4b: Single-tone ICI leakage ---
k0 = 10;   % data-subcarrier index

Xmod_tone = zeros(1, Nused);
Xmod_tone(k0) = 1;
X_shift_tone = map_to_subcarriers(Xmod_tone, Nfft, Nused, Nvc);

% One representative channel realization
h = zeros(1, L);
for l = 1:L
    h(l) = sqrt(pdp(l)) * Rayleigh_MIMO(1,1);
end

x_cp   = ofdm_mod(X_shift_tone, Ng);
x_core = ifft(X_shift_tone);
x_zp   = [x_core, zeros(1, Ng)];

y_cp = conv(x_cp, h);
y_cp = y_cp(1:Nfft+Ng);

y_zp = conv(x_zp, h);
y_zp = y_zp(1:Nfft+Ng);

Y_cp  = fft(y_cp(Ng+1 : Ng+Nfft));
Y_zpn = fft(y_zp(1:Nfft));

r_oa = y_zp(1:Nfft);
r_oa(1:Ng) = r_oa(1:Ng) + y_zp(Nfft+1 : Nfft+Ng);
Y_zpo = fft(r_oa);

figure('Name', 'Section 4: CP vs ZP', 'Position', [100 100 1200 500]);

subplot(1,2,1);
semilogy(1:Nused, mse_cp,  'b-o', 'LineWidth', 1.4, 'MarkerSize', 4); hold on;
semilogy(1:Nused, mse_zpn, 'r-s', 'LineWidth', 1.4, 'MarkerSize', 4);
semilogy(1:Nused, mse_zpo, 'g-^', 'LineWidth', 1.4, 'MarkerSize', 4);
grid on;
xlabel('Data subcarrier');
ylabel('MSE');
title(sprintf('Per-subcarrier MSE, noiseless multipath (L=%d, Ng=%d)', L, Ng));
legend('CP', 'ZP-naive', 'ZP+OA', 'Location', 'best');

subplot(1,2,2);
stem(0:Nfft-1, abs(Y_cp),  'b-o', 'LineWidth', 1.2); hold on;
stem(0:Nfft-1, abs(Y_zpn), 'r-s', 'LineWidth', 1.2);
stem(0:Nfft-1, abs(Y_zpo), 'g-^', 'LineWidth', 1.2);
grid on;
xlabel('FFT bin k');
ylabel('|Y[k]|');
title(sprintf('Single-tone Rx spectrum (tone on data bin %d)', k0));
legend('CP', 'ZP-naive', 'ZP+OA', 'Location', 'best');

%% --- Local function ----------------------------------
function Hdat = H_data_bins(H, Nfft, Nused, Nvc)
% Pick the data-subcarrier entries of a length-Nfft frequency response.
    Hdat = [H(Nfft-Nused/2+1:Nfft), H(2:Nused/2+1)];
end


