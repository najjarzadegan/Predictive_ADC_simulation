clc

step = 3

p0_1 = -0.7;
p1_1 = 0.9;
thp_1 = 105;
z0_1 = 0.9;
z1_1 = 0.92;
thz_1 = 35;

if(step == 1)

p0 = -0.9:0.3:0.9;
p1 = 0.05:0.1:0.95;
thp = 0:45:180;

z0 = -0.9:0.3:0.9;
z1 = 0.05:0.1:0.95;
thz = 0:45:180;

elseif(step == 2)

p0 = max(p0_1-0.2,-0.99):0.05:min(p0_1+0.2,0.99);
p1 = max(p1_1-0.07,-0.99):0.02:min(p1_1+0.07,0.99);
thp = thp_1-30:10:thp_1+30;

z0 = max(z0_1-0.2,-0.99):0.05:min(z0_1+0.2,0.99);
z1 = max(z1_1-0.07,-0.99):0.02:min(z1_1+0.07,0.99);
thz = thz_1-30:10:thz_1+30;

elseif(step == 3)

p0 = p0_1;
p1 = max(p1_1-0.02,-0.99):0.01:min(p1_1+0.02,0.99);
thp = thp_1-7:1:thp_1+7;

z0 = z0_1;
z1 = max(z1_1-0.02,-0.99):0.01:min(z1_1+0.02,0.99);
thz = thz_1-7:1:thz_1+7;

end

DAC_r = 10
ADC_r = 10
K = 11
OSR = 5

count = 0;

Emax = zeros([length(p0),length(p1),length(thp),length(z0),length(z1),length(thz)]);
prod(size(Emax))


for i1 = 1:length(p0)
    for i2 = 1:length(p1)
        for i3 = 1:length(thp)
            for i4 = 1:length(z0)
                for i5 = 1:length(z1)
                    for i6 = 1:length(thz)

                        [Emax(i1,i2,i3,i4,i5,i6),~] = CLF_optimizer(DAC_r,ADC_r,K,OSR,z0(i4),z1(i5),thz(i6),p0(i1),p1(i2),thp(i3));
                        count = count+1;

                        if(rem(count,10000) == 0)
                            count
                        end
                    end
                end
            end
        end
    end
end

size(Emax)
[m,I] = min(Emax,[],"all")

indx6 = floor((I-1)/(length(p0)*length(p1)*length(thp)*length(z0)*length(z1)))+1
I = I-(indx6-1)*(length(p0)*length(p1)*length(thp)*length(z0)*length(z1))

indx5 = floor((I-1)/(length(p0)*length(p1)*length(thp)*length(z0)))+1
I = I-(indx5-1)*(length(p0)*length(p1)*length(thp)*length(z0))

indx4 = floor((I-1)/(length(p0)*length(p1)*length(thp)))+1
I = I-(indx4-1)*(length(p0)*length(p1)*length(thp))

indx3 = floor((I-1)/(length(p0)*length(p1)))+1
I = I-(indx3-1)*(length(p0)*length(p1))

indx2 = floor((I-1)/(length(p0)))+1
I = I-(indx2-1)*(length(p0))

indx1 = I


Emax(indx1,indx2,indx3,indx4,indx5,indx6)

param = [p0(indx1),p1(indx2),thp(indx3),z0(indx4),z1(indx5),thz(indx6)]

[a, b] = CLF_optimizer(DAC_r,ADC_r,K,OSR,z0(indx4),z1(indx5),thz(indx6),p0(indx1),p1(indx2),thp(indx3))

zero = [z0(indx4) z1(indx5)*(cos(pi*thz(indx6)/180)+1i*sin(pi*thz(indx6)/180)) z1(indx5)*(cos(pi*thz(indx6)/180)-1i*sin(pi*thz(indx6)/180))];
pole = [p0(indx1) p1(indx2)*(cos(pi*thp(indx3)/180)+1i*sin(pi*thp(indx3)/180)) p1(indx2)*(cos(pi*thp(indx3)/180)-1i*sin(pi*thp(indx3)/180))];

close all
title1 = sprintf('M = %d, N = %d, K = %d, OSR = %d',DAC_r,ADC_r,K,OSR);
zplane(zero',pole')
title(title1)

figure
Predictor_Coeff(0, 0, 0, 0, 0, 0, K,0,pole,zero)
