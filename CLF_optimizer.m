function [Emax,portion] = CLF_optimizer(DAC_r,ADC_r,K,OSR,z0,z1,thz,p0,p1,thp)


Ap = 1;
theta = pi/OSR;

zero = [z0 z1*(cos(pi*thz/180)+1i*sin(pi*thz/180)) z1*(cos(pi*thz/180)-1i*sin(pi*thz/180))];
pole = [p0 p1*(cos(pi*thp/180)+1i*sin(pi*thp/180)) p1*(cos(pi*thp/180)-1i*sin(pi*thp/180))];

%%%%%%%% Filter Generation %%%%%%%%%%%%%%%%%%%%%%
L1 = 0;
L2 = 0;
 
[N,D] = zp2tf(zero',pole',1);
       
H_N1 = (D - N)/K
H_D1 = N
 
    N_noise = K*conv(N,H_N1);
    D_noise = conv(D,H_D1);
    
    try
        L2 = filternorm(N_noise,D_noise);
    catch ME
        warning('The filter is almost unstable')
    end

%     try
%         L22 = filternorm(D-N,D);
%     catch ME
%         warning('The filter is almost unstable')
%     end
    
    [h_loop,w_loop] = freqz(N,D);
    [h_H,w_H] = freqz(H_N1,H_D1);
    
    
    L1 = max(abs(h_loop(1:floor(length(w_loop)/OSR))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Term2 = L1;
Term3 = L2;
%L22;
M0 = log2((1 + L2)./(1./K' - 1*L1));
DoC = (1./K')./(1./K' - 1*L1); % Degree of Confidence
Term = [L1,L2,M0,DoC];

Emax = 2^(-DAC_r) + Ap*L1 + L2*sqrt(3*(4^(-ADC_r)+K^2*4^(-DAC_r)))/K;
portion = [2^(-DAC_r) Ap*L1 L2*sqrt(3*(4^(-ADC_r)+K^2*4^(-DAC_r)))/K]*(100/Emax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title1 = sprintf('M = %d, N = %d, K = %d, OSR = %d',DAC_r,ADC_r,K,OSR);
zplane(zero',pole')
title(title1)

end