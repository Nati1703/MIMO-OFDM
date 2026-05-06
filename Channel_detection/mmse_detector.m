function x_hat = mmse_detector(y, H, noise_var)
% MMSE MIMO detection
%
% Inputs:
%   y         - received signal (Nr x 1)
%   H         - channel matrix (Nr x Nt)
%   noise_var - noise variance sigma^2
%
% Output:
%   x_hat - detected signal (Nt x 1), continuous-valued (soft estimate)
%
% Formula: x_hat = (H^H*H + sigma^2*I)^{-1} * H^H * y

    Nt = size(H, 2);
    x_hat = (H' * H + noise_var * eye(Nt)) \ (H' * y);
end
