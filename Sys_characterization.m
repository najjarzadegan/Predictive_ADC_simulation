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

Out(:,i) = Calib(:,3*(i-1)+1);
In(:,i)  = Calib(:,3*(i-1)+3);
Noise(:,i) = Calib(:,3*(i-1)+2);

end

C_Noise2 = sum(Noise(:,1).^2)/l
C_Noise3 = sum(Noise(:,1).^3)/l
C_Noise4 = sum(Noise(:,1).^4)/l

C_Noise2 = sum(sum(Noise.^2))/(l_a*l)
C_Noise3 = sum(sum(Noise.^3))/(l_a*l)
C_Noise4 = sum(sum(Noise.^4))/(l_a*l)


C_Out_Noise = zeros(l_a,1);
C_In_Noise2 = zeros(l_a,1);
C_In2_Noise2 = zeros(l_a,1);
C_In_Noise3 = zeros(l_a,1);

C_In_N = zeros(l_a,1);
C_In2_N = zeros(l_a,1);
C_In3_N = zeros(l_a,1);

for i = 1:l_a

C_Out_Noise(i) = sum(Out(:,i).*Noise(:,i))/l;
C_In_Noise2(i) = sum(In(:,i).*Noise(:,i).^2)/l;
C_In2_Noise2(i) = sum(In(:,i).^2.*Noise(:,i).^2)/l;
C_In_Noise3(i) = sum(In(:,i).*Noise(:,i).^3)/l;

C_In_N(i) = sum(In(:,i).*Noise(:,i))/l
C_In2_N(i) = sum(In(:,i).^2.*Noise(:,i))/l
C_In3_N(i) = sum(In(:,i).^3.*Noise(:,i))/l

end

% C_Out_Noise(5x1) = Matrix(5x3)*C(3x1)
% Matrix = [C_Noise2*ones(5,1) C_Noise3*ones(5,1)+2*C_In_Noise2
% C_Noise4*ones(5,1)+3*C_In_Noise3+3*C_In2_Noise2]
A1 = 1e3;
A2 = 1e-1;
A3 = 1e-4;
A4 = 1e-8;

Matrix = [A1*mean(mean(Noise))*ones(l_a,1) A2*(C_Noise2*ones(l_a,1)+C_In_N) A3*(C_Noise3*ones(l_a,1)+2*C_In_Noise2+C_In2_N) A4*(C_Noise4*ones(l_a,1)+3*C_In_Noise3+3*C_In2_Noise2+C_In3_N)]

C_e = (Matrix'*Matrix)\(Matrix'*C_Out_Noise);

Ce = C_e'.*[A1 A2 A3 A4]

Ce_inv = inv([1 Ce(1) Ce(1)^2 Ce(1)^3;0 Ce(2) 2*Ce(1)*Ce(2) 3*Ce(1)*Ce(2);0 Ce(3) 2*Ce(1)*Ce(3)+Ce(2)^2 3*Ce(1)*(Ce(1)*Ce(3)+Ce(2)^2);0 Ce(4) 2*Ce(1)*Ce(4)+2*Ce(2)*Ce(3) 3*Ce(1)^2*Ce(4)+6*Ce(1)*Ce(2)*Ce(3)+Ce(2)^3])*[0;1;0;0]
