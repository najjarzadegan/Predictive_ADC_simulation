function [N,D] = Predictor_Coeff(type, alpha, beta, order, angle, OSR, gain,normalized,P,Z)

%close all
z = tf('z');

n = floor(order/2)
teta0 = pi/angle/OSR/n

if(rem(order,2) == 0)
    vec1 = 1:n;
    vec2 = -vec1;
    vec = [flip(vec2) vec1]
else
    vec = -n:n
end

teta  = vec*teta0
zero  = alpha*exp(1i*teta)
pole  = beta*exp(1i*teta)

if(type == 1)

    if(normalized == 1)
        k   = (beta+1)^order/(alpha+1)^order
    else
        k = 1
    end

    clf = zpk(alpha*ones(1,order),beta*ones(1,order),k,1)

elseif(type == 2)

    if(normalized == 1)
        k   = (beta+1)^order/abs(prod(1+zero))
    else
        k = 1
    end

    clf = zpk(zero,beta*ones(1,order),k,1)

elseif(type == 3)

    if(normalized == 1)
        k   = abs(prod(1+pole))/abs(prod(1+zero))
    else
        k = 1
    end

    clf = zpk(zero,pole,k,1)

elseif(type == 0)

    clf = zpk(Z,P,1,1)

else
    disp('Filter type is not correct.')
end

H=(1/clf-1)*z/gain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_h, den_h] = tfdata(H, 'v')

Y_coeff = -den_h(2:end)/den_h(1)
X_coeff = num_h/den_h(1)

[Mag_h,w_h] = freqz(num_h, den_h,1000);
plot(w_h,abs(Mag_h))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num_clf, den_clf] = tfdata(clf, 'v');

[Mag_clf,w_clf] = freqz(num_clf, den_clf,1000);
figure
plot(w_clf,abs(Mag_clf))


