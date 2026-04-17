function H = Rician_model(K_dB,L)
    % Rician channel with LoS phase
    %Input :
    %   K_dB : Rician K-factor in dB
    %   L:  number of samples
    K = 10 ^(K_dB/10);
%     phi = 2 * pi * rand; %random phase as well% lets put it to zero for now
%     phi = 0;
    H = sqrt(K/(K+1))+ sqrt(1/(K+1)) * Rayleigh_model(L);

end
