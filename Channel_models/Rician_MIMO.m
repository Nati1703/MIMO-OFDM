function H = Rician_MIMO(Nr,Nt,K_dB)
H = zeros(Nr * Nt);
for rx = 1:Nr
    for tx = 1:Nt
        H(rx,tx) = Rician_model(K_dB,1);
    end 
end

end