function Hls = ls_ce(Y, X)
%LS_CE Least-squares channel estimate on known tones.
    Xsafe = X;
    Xsafe(abs(Xsafe) < 1e-12) = 1e-12;

    Hls = Y ./ Xsafe; % see equation 6.6 from cho's book.
    % Simplified for where X is just per subcarrier pilot right.
    % General LS in NLASP assignmets (X^HX)^-1X^HY. that woudl be the MIMO
    % equivalent; also like if you can't inverse you know the stuff, pinv
    % (See Ivlrac's notes)
end
