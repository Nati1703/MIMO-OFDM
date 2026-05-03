function Xhat = one_tap_eq(Y, H)
%ZF equalizer. 
%ONE_TAP_EQ One-tap frequency-domain equalizer.
%
%   Xhat = one_tap_eq(Y, H)

    if ~isequal(size(Y), size(H))
        error('Y and H must have the same size.');
    end

    Hsafe = H;
    Hsafe(abs(Hsafe) < 1e-10) = 1e-10;

    Xhat = Y ./ Hsafe;
end