function Hmmse = mmse_ce(Hls, noise_var, pdp, Nfft, Nvc)
%MMSE_CE LMMSE denoising of a frequency-domain LS channel estimate.
%
%   Hmmse = mmse_ce(Hls, noise_var, pdp)
%       Legacy. Builds R assuming the input sits on consecutive bins
%       0..Nused-1 with Nfft_corr = Nused. Biased when there are virtual
%       carriers (the real bins are not consecutive). Kept for back-compat.
%
%   Hmmse = mmse_ce(Hls, noise_var, pdp, Nfft_corr)
%       Legacy with a user-supplied Nfft_corr scalar. Still wrong delta.
%
%   Hmmse = mmse_ce(Hls, noise_var, pdp, Nfft, Nvc)
%       Recommended. Recovers the actual used-subcarrier bins k_m and
%       builds
%             R(m1, m2) = sum_l pdp(l) * exp(-j*2*pi*l*(k_m1 - k_m2)/Nfft)
%       so the correlation model matches the physical spectrum. Then
%             Hmmse = R * (R + noise_var * I)^-1 * Hls.
%       At high SNR this asymptotes to the L-tap projector (matches
%       DFT-LS); at low SNR it adds noise-aware shrinkage on top.

    was_row = isrow(Hls);
    hls = Hls(:);
    Nused = numel(hls);

    if nargin < 3 || isempty(pdp)
        pdp = 1;
    end
    pdp = pdp(:);
    pdp = pdp / sum(pdp);

    full_mode = (nargin >= 5) && ~isempty(Nfft) && ~isempty(Nvc);

    if full_mode
        % Recover actual FFT-bin index of each input position by probing
        % map_to_subcarriers with a uniquely-valued vector.
        probe = 1:Nused;
        probe_full = map_to_subcarriers(probe, Nfft, Nused, Nvc);
        probe_int = round(real(probe_full));
        k_used = zeros(Nused, 1);
        for m = 1:Nused
            pos = find(probe_int == m, 1);
            k_used(m) = pos - 1;        % 0-indexed bin
        end
        Nfft_corr = Nfft;
    else
        % Legacy: assume consecutive bins, optional Nfft_corr in arg 4.
        k_used = (0:Nused-1).';
        if nargin >= 4 && ~isempty(Nfft)
            Nfft_corr = Nfft;
        else
            Nfft_corr = Nused;
        end
    end

    delta = k_used - k_used.';

    R = zeros(Nused, Nused);
    for ell = 0:numel(pdp)-1
        R = R + pdp(ell+1) * exp(-1j * 2*pi * ell * delta / Nfft_corr);
    end

    Hmmse = R / (R + noise_var * eye(Nused)) * hls;

    if was_row
        Hmmse = Hmmse.';
    end
end
