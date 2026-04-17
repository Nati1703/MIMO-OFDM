function H = Rayleigh_model(L)
% Generates Rayleigh channel coefficients
% Number of channel realizations
H = ((randn(1,L) + j* randn(1,L))/sqrt(2));
end