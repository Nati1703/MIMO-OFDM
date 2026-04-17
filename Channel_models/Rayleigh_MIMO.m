function H = Rayleigh_MIMO(Nr,Nt)
H = Rayleigh_model(Nr*Nt);
H = reshape(H, Nr, Nt);
end
