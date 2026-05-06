function x_hat = zf_detector(y, H)
% Zero-forcing MIMO detection
%
% Inputs:
%   y - received signal (Nr x 1)
%   H - channel matrix (Nr x Nt)
%
% Output:
%   x_hat - detected signal (Nt x 1), continuous-valued (soft estimate)
%
% Formula: x_hat = (H^H * H)^{-1} * H^H * y = H^+ * y
% Same as LS estimation — just inverts the channel
% Problem: amplifies noise when H is ill-conditioned

    x_hat = (H' * H) \ (H' * y);
end
