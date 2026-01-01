clc

ADC1 = readtable('temp6_1kHz500mV_7bitDAC_8bitADC_Noise_1.txt');
ADC2 = readtable('temp6_1kHz600mV_7bitDAC_8bitADC_Noise_1.txt');
ADC3 = readtable('temp6_1kHz800mV_7bitDAC_8bitADC_Noise_1.txt');
ADC4 = readtable('temp6_1kHz1000mV_7bitDAC_8bitADC_Noise_1.txt');
ADC5 = readtable('temp6_1kHz1100mV_7bitDAC_8bitADC_Noise_1.txt');
ADC6 = readtable('temp6_1kHz1200mV_7bitDAC_8bitADC_Noise_1.txt');
ADC7 = readtable('temp6_1kHz1300mV_7bitDAC_8bitADC_Noise_1.txt');

A = [500,600,800,1000,1100,1200,1300];

l_a = length(A);


l = height(ADC1);
Xr = zeros(1,l);

Calib = zeros(l,21);
Calib(:,1:3) = table2array(ADC1);
Calib(:,4:6) = table2array(ADC2);
Calib(:,7:9) = table2array(ADC3);
Calib(:,10:12) = table2array(ADC4);
Calib(:,13:15) = table2array(ADC5);
Calib(:,16:18) = table2array(ADC6);
Calib(:,19:21) = table2array(ADC7);

for i = 1:l_a

ADC(:,i) = Calib(2:4096,3*(i-1)+1);
DAC(:,i)  = Calib(1:4095,3*(i-1)+3) + Calib(1:4095,3*(i-1)+2);

end

ADCa = ADC(:,7);
DACa = DAC(:,7);

C = [2300, 11.56, 0, 0];

f_DAC = C(1) + C(2)*DACa + C(3)*DACa.^2 + C(4)*DACa.^3;


[p_inv1, out1, p_inv2, out2, out2_check] = Poly_inv(C,f_DAC + ADCa);

norm(out2-out2_check)

ENOB = (sinad(round(out2),50e3)-1.76)/6.02