clc
clear all

fin = 100;
T = 1e-4;
P = 1e5;
time = 0:T:(P-1)*T;


I0 = 5;
I1 = 1;
I2 = 0.0001;

I = I0 + I1*sin(2*pi*fin*time) + I2*sin(4*pi*fin*time);

H0 = sum(I)/P

H1a = 2*sum(I.*sin(2*pi*fin*time))/P;
H1b = 2*sum(I.*cos(2*pi*fin*time))/P;
H1  = sqrt(H1a^2 + H1b^2)

H2a = 2*sum(I.*sin(4*pi*fin*time))/P;
H2b = 2*sum(I.*cos(4*pi*fin*time))/P;
H2  = sqrt(H2a^2 + H2b^2)

H3a = 2*sum(I.*sin(6*pi*fin*time))/P;
H3b = 2*sum(I.*cos(6*pi*fin*time))/P;
H3  = sqrt(H3a^2 + H3b^2)

THD2 = 20*log10(H1/H2)




        H2a_INA_1 = p0*Ts/(1+p0*Ts)*2*sin(4*pi*fin*time(i))*Vs_a(sample) + H2a_INA_1/(1+p0*Ts);
        H2a_INA_2 = p0*Ts/(1+p0*Ts)*H2a_INA_1 + H2a_INA_2/(1+p0*Ts);
        H2a_INA_3 = p0*Ts/(1+p0*Ts)*H2a_INA_2 + H2a_INA_3/(1+p0*Ts);
        H2a_vec(sample) = H2a_INA_3;
        
        H2b_INA_1 = p0*Ts/(1+p0*Ts)*2*cos(4*pi*fin*time(i))*Vs_a(sample) + H2b_INA_1/(1+p0*Ts);
        H2b_INA_2 = p0*Ts/(1+p0*Ts)*H2b_INA_1 + H2b_INA_2/(1+p0*Ts);
        H2b_INA_3 = p0*Ts/(1+p0*Ts)*H2b_INA_2 + H2b_INA_3/(1+p0*Ts);
        H2b_vec(sample) = H2b_INA_3;
        
        H2_INA(sample)  = sqrt(H2a_INA_3^2 + H2b_INA_3^2);