function Hhat = dft_ce(Hin, L, Nfft, Nvc)
%DFT_CE Project a noisy frequency-domain channel estimate onto the
%       L-tap channel subspace.
%
%   Hhat = dft_ce(Hin, L)
%       Legacy mode. Truncates in the numel(Hin)-point delay domain,
%       which biases the result whenever Hin is the used-subcarrier
%       slice of a larger Nfft spectrum (because the channel taps are
%       sparse only in the Nfft-point IFFT). Kept for back-compat.
%
%   Hhat = dft_ce(Hin, L, Nfft, Nvc)
%       Recommended. Solves
%             h_est = argmin_h || Hin - A*h ||^2
%       with A(m, l+1) = exp(-j*2*pi*k_m*l/Nfft) on the *actual*
%       used-subcarrier indices k_m, then returns A*h_est. This is the
%       orthogonal projector onto the L-tap channel subspace, so it has
%       no spectral-leakage floor. The expected MSE is L/Nused times
%       the LS MSE.
%
% Used-subcarrier indices are recovered by probing map_to_subcarriers
% with a uniquely-valued vector, so this works for any mapping
% convention that map_to_subcarriers / demap_from_subcarriers implement.

    was_row = isrow(Hin);
    if was_row
        Hrow = Hin;
    else
        Hrow = Hin.';
    end
    Nused = numel(Hrow);

    full_mode = (nargin >= 4) && ~isempty(Nfft) && ~isempty(Nvc);

    if full_mode
        % Recover the actual FFT-bin index of each input position by
        % sending a unique-valued probe through map_to_subcarriers.
        probe = 1:Nused;
        probe_full = map_to_subcarriers(probe, Nfft, Nused, Nvc);
        probe_int = round(real(probe_full));        % robust to tiny FP noise
        k_used = zeros(1, Nused);
        for m = 1:Nused
            pos = find(probe_int == m, 1);
            k_used(m) = pos - 1;                    % convert to 0-indexed bin
        end

        % L-tap DFT basis on the used bins
        l_idx = 0:L-1;
        A = exp(-1j * 2*pi * (k_used(:) * l_idx) / Nfft);   % Nused x L

        % LS estimate of channel taps, then project back
        h_est = A \ Hrow(:);
        Hhat_row = (A * h_est).';
    else
        if L >= Nused
            Hhat_row = Hrow;
        else
            h = ifft(Hrow);
            h(L+1:end) = 0;
            Hhat_row = fft(h);
        end
    end

    if was_row
        Hhat = Hhat_row;
    else
        Hhat = Hhat_row.';
    end
end
