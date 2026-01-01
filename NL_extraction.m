% nonlinear coefficient calculation

clc
close all
f = 1;
OSR = 10;
N = 4000;
C = [0,2,0.01,0.001];

dt = 0.5/f/OSR;
T = N/f;
t = 0:dt:T;
l = length(t);

nn = generate_noise(3, 0, l);
Noise = nn';

C_Noise2 = sum(Noise.^2)/l
C_Noise3 = sum(Noise.^3)/l
C_Noise4 = sum(Noise.^4)/l


A = [50,60,70,80,90,100];
%A = [80,90,100];
l_a = length(A);

C_Out_Noise = zeros(l_a,1);
C_In_Noise2 = zeros(l_a,1);
C_In2_Noise2 = zeros(l_a,1);
C_In_Noise3 = zeros(l_a,1);

C_In_N = zeros(l_a,1);
C_In2_N = zeros(l_a,1);
C_In3_N = zeros(l_a,1);

for i = 1:l_a

In = A(i)*sin(2*pi*f*t);
X = In + Noise;

% assume Y = C0 + C1*X + C2*X^2 + C3*X^3 and X = X0 + n 
% So, Y = C0 + C1*X0 + C2*X0^2 + C3*X0^3
%       + n*(C1 + 2*C2*X0 + 3*C3*X0^2)
%       + n^2*(C2 + 3*C3*X0)
%       + n^3*(C3)

Out = C(1) + C(2)*X + C(3)*X.^2 + C(4)*X.^3;

% <Y,n> = 0 + 0 + 0 + 0
%       + C1*<n,n> + 2*C2*<X0,n^2> + 3*C3*<X0^2,n^2>
%       + C2*<n,n^2> + 3*C3*<X0,n^3>
%       + C3*<n,n^3>

% C_A_B = <A,B> = sum(A.*B)/length(A)

C_Out_Noise(i) = sum(Out.*Noise)/l;
C_In_Noise2(i) = sum(In.*Noise.^2)/l;
C_In2_Noise2(i) = sum(In.^2.*Noise.^2)/l;
C_In_Noise3(i) = sum(In.*Noise.^3)/l;

C_In_N(i) = sum(In.*Noise)/l
C_In2_N(i) = sum(In.^2.*Noise)/l
C_In3_N(i) = sum(In.^3.*Noise)/l

end

% C_Out_Noise(5x1) = Matrix(5x3)*C(3x1)
% Matrix = [C_Noise2*ones(5,1) C_Noise3*ones(5,1)+2*C_In_Noise2
% C_Noise4*ones(5,1)+3*C_In_Noise3+3*C_In2_Noise2]
A1 = 1e4;
A2 = 1e1;
A3 = 1;
A4 = 1e-3;

Matrix = [A1*mean(Noise)*ones(l_a,1) A2*(C_Noise2*ones(l_a,1)+C_In_N) A3*(C_Noise3*ones(l_a,1)+2*C_In_Noise2+C_In2_N) A4*(C_Noise4*ones(l_a,1)+3*C_In_Noise3+3*C_In2_Noise2+C_In3_N)]

C_e = (Matrix'*Matrix)\(Matrix'*C_Out_Noise);

Ce = C_e'.*[A1 A2 A3 A4]

Ce_inv = inv([1 Ce(1) Ce(1)^2 Ce(1)^3;0 Ce(2) 2*Ce(1)*Ce(2) 3*Ce(1)*Ce(2);0 Ce(3) 2*Ce(1)*Ce(3)+Ce(2)^2 3*Ce(1)*(Ce(1)*Ce(3)+Ce(2)^2);0 Ce(4) 2*Ce(1)*Ce(4)+2*Ce(2)*Ce(3) 3*Ce(1)^2*Ce(4)+6*Ce(1)*Ce(2)*Ce(3)+Ce(2)^3])*[0;1;0;0]

semilogy(abs(fft(Out)))
figure
semilogy(abs(fft(In)))


