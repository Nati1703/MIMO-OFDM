function x_hat = osic_detector(y, H, noise_var, constellation)
% Ordered Successive Interference Cancellation (OSIC) detection
% Uses MMSE weights with SINR-based ordering
%
% Inputs:
%   y             - received signal (Nr x 1)
%   H             - channel matrix (Nr x Nt)
%   noise_var     - noise variance sigma^2
%   constellation - vector of valid constellation points (e.g., qammod(0:M-1,M,...))
%
% Output:
%   x_hat - detected signal (Nt x 1), hard decisions
%
% Algorithm:
%   1. Compute MMSE weights for all remaining streams
%   2. Pick stream with highest post-detection SINR
%   3. Detect it (hard decision → nearest constellation point)
%   4. Subtract its contribution from received signal
%   5. Repeat with reduced system

    [Nr, Nt] = size(H);
    x_hat = zeros(Nt, 1);
    remaining = 1:Nt;
    y_res = y;
    H_res = H;

    for stage = 1:Nt
        Nt_rem = length(remaining);

        % MMSE weight matrix for remaining streams
        W = (H_res' * H_res + noise_var * eye(Nt_rem)) \ H_res';

        % Compute post-detection SINR for each remaining stream
        sinr = zeros(Nt_rem, 1);
        for i = 1:Nt_rem
            wi = W(i, :).';
            % Signal power: |w_i^H h_i|^2
            sig = abs(wi' * H_res(:, i))^2;
            % Interference: sum over other streams
            interf = 0;
            for j = [1:i-1, i+1:Nt_rem]
                interf = interf + abs(wi' * H_res(:, j))^2;
            end
            % Noise: sigma^2 * ||w_i||^2
            noi = noise_var * (wi' * wi);
            sinr(i) = sig / (interf + noi);
        end

        % Pick the stream with highest SINR
        [~, best_idx] = max(sinr);
        best_stream = remaining(best_idx);

        % Detect: MMSE estimate → snap to nearest constellation point
        x_soft = W(best_idx, :) * y_res;
        x_hat(best_stream) = slice_to_constellation(x_soft, constellation);

        % Cancel detected stream from received signal
        y_res = y_res - H_res(:, best_idx) * x_hat(best_stream);

        % Remove detected column from H and tracking
        H_res(:, best_idx) = [];
        remaining(best_idx) = [];
    end
end

function x_nearest = slice_to_constellation(x_soft, constellation)
% Snap complex value to nearest constellation point
    [~, idx] = min(abs(x_soft - constellation));
    x_nearest = constellation(idx);
end
