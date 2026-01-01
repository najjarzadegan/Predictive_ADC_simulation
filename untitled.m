clc

% zero = [z0(indx4) z1(indx5)*(cos(pi*thz(indx6)/180)+1i*sin(pi*thz(indx6)/180)) z1(indx5)*(cos(pi*thz(indx6)/180)-1i*sin(pi*thz(indx6)/180))];
% pole = [p0(indx1) p1(indx2)*(cos(pi*thp(indx3)/180)+1i*sin(pi*thp(indx3)/180)) p1(indx2)*(cos(pi*thp(indx3)/180)-1i*sin(pi*thp(indx3)/180))];

z0 = 0.8
z1 = 0.92
thz = 33

p0=-0.4
p1=0.9
thp=69

zero = [z0 z1*(cos(pi*thz/180)+1i*sin(pi*thz/180)) z1*(cos(pi*thz/180)-1i*sin(pi*thz/180))];
pole = [p0 p1*(cos(pi*thp/180)+1i*sin(pi*thp/180)) p1*(cos(pi*thp/180)-1i*sin(pi*thp/180))];


% [CLF_n,CLF_d]=zp2tf(zero',pole',1)
% 
% sys = tf(CLF_n,CLF_d);
% pzmap(sys)

zplane(zero',pole')

hold on

z0 = 0.99
z1 = 0.99
thz = 36

p0=0.8
p1=0.8
thp=36

zero = [z0 z1*(cos(pi*thz/180)+1i*sin(pi*thz/180)) z1*(cos(pi*thz/180)-1i*sin(pi*thz/180))];
pole = [p0 p1*(cos(pi*thp/180)+1i*sin(pi*thp/180)) p1*(cos(pi*thp/180)-1i*sin(pi*thp/180))];

zplane(zero',pole','b')







DAC_r = 6
ADC_r = 6
K = 11.56
OSR = 5

%Predictor_Coeff(0, 0, 0, 0, 0, 0, K,0,pole,zero)

%[a, b] = CLF_optimizer(DAC_r,ADC_r,K,OSR,z0,z1,thz,p0,p1,thp)
