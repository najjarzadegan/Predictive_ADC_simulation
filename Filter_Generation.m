% Filter generation
clc
close all

G = 1
Order = 3
theta = pi/16
K = 16
a = 0.95
b = 0
Type = 2

if Type == 1
%%%%%%  Type = 1  %%%%%%%%%

p = b*ones(Order,1);
z = a*ones(Order,1);

elseif Type == 2
%%%%%%  Type = 2  %%%%%%%%%

z = 0*ones(Order,1);
p = 0*ones(Order,1);

if mod(Order,2) == 1
    z(Order) = a; 
end

NZ = floor(Order/2);
if NZ > 0
for i = 1:NZ
    
    A = roots([1 -2*a*cos(theta/NZ) a^2]);
    z(2*i-1) = A(1);
    z(2*i) = A(2);
end
end

elseif Type == 3
%%%%%%  Type = 3  %%%%%%%%%

z = 0*ones(Order,1);
p = 0*ones(Order,1);

if mod(Order,2) == 1
    z(Order) = a;
    p(Order) = b;
end

NZ = floor(Order/2);
if NZ > 0
for i = 1:NZ
    
    A = roots([1 -2*a*cos(theta/NZ) a^2]);
    z(2*i-1) = A(1);
    z(2*i) = A(2);
    
    A = roots([1 -2*b*cos(theta/NZ) b^2]);
    p(2*i-1) = A(1);
    p(2*i) = A(2);
    
end 
end
    
end    
    
[N,D] = zp2tf(z,p,G);

H_N = (D - N)/K;
H_D = N;

In = ADC;
Out = 0*ADC;

for i = 1:length(In)
    if i < Order + 1
        Out(i) = 0;
    else
%         Out(i - Order:i - 1)
%         (flip(-H_D(2:end)))'
%         In(i - Order + 1:i)
%         (flip(H_N(2:end)))'
    Out(i) = Out(i - Order:i - 1)*(flip(-H_D(2:end)))' + ...
             In(i - Order + 1:i)*(flip(H_N(2:end)))';
    end
end

    