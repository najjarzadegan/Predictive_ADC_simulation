ADC = readtable('temp10_DC1100mVoffset0V_KEITHLEY2602_7bitDAC_9bitADC.txt');
l = height(ADC);
Xr = zeros(1,l);
ADC = table2array(ADC);

%gain = 11.65;
%gain = 22.29;
gain = 28.44;

fs   = 50e3;
k = -0;


freq = 0:fs/l:(l-1)*fs/l;
close all

adc = ADC(2:l,1);
dac = ADC(1:l-1,3);

R = 4.82e3;
C = 3.3e-9;
T = 20e-6;

alpha = R*C/T;
adc_m = (1+alpha)*adc - alpha*[0;adc(1:l-2)];

for i = 2:l

    Xr12(i) = ADC(i,1)/gain + ADC(i-1,2);
    Xr7(i) = ADC(i,1)/gain + ADC(i-1,3) - k*ADC(i-1,2);
    %Xr7(i) = adc_m(i-1)/gain + ADC(i-1,3);

end


%Xr7 = adc_m/gain + dac;

% plot(Xr7)
% title('Reconstructed input')
% xlabel('Sample')
% ylabel('Digital Code')
% ylim([0,4095])
% 
% figure
% plot(ADC(:,1))
% title('ADC')
% xlabel('Sample')
% ylabel('Digital Code')
% ylim([0,4095])
% 
% figure
% plot(ADC(:,3))
% title('ADC')
% xlabel('Sample')
% ylabel('Digital Code')
% ylim([0,4095])
% 
% figure
% F = abs(fft(Xr7)/4096);
% plot(freq(1:l/2),20*log10(F(1:l/2)))
% title('Input FFT')
% xlabel('Frequency')
% ylabel('dB')
% 
% figure
% 
% amp = 1.3
% Afs = 1.5;
% 
% 
% 
% xm = 2:4096;
% %gs=outputparameters(1)*sin(outputparameters(2)*xm+outputparameters(3))+outputparameters(4);
% 
% 
% sinad(Xr7,fs)
% 
% figure
% sinad(adc,fs)
% title('adc')
% figure
% sinad(dac,fs)
% title('dac')
% 
% (sinad(Xr7,fs)-1.76)/6.02
% (sinad(Xr7,fs)-1.76)/6.02+log2(Afs/amp)

Avg = mean(Xr7)

%log2(gain)
% figure
% plot(gs,dac)
