function [FF_N,FF_D,IF_N,IF_D,freq_FF,Resp_FF,freq_IF,Resp_IF] = ...
    FF_IF_Designer(FF_Order,FF_Limit,Est_FF,Est_IF,FF_Type1,N,OSR)
clc
close all
N = str2num(N);
FF_N = {};
FF_D = {};
IF_N = {};
IF_D = {};
freq_FF = {};
Resp_FF = {};
freq_IF = {};
Resp_IF = {};
for i = 1:N-1
    
    if ~Est_FF(2,i) && Est_FF(1,i)
        
        if strcmp(FF_Type1,'Butterworth')
            
            [FF_N{i},FF_D{i}] = butter(FF_Order(i),1/OSR,'low');
            
        elseif strcmp(FF_Type1,'FIR')
            
            [FF_N{i},FF_D{i}] = fir1(FF_Order(i),1/OSR);
            
        elseif strcmp(FF_Type1,'Chebyshev I')
            
            [FF_N{i},FF_D{i}] = cheby1(FF_Order(i),10,1/OSR);   % 10 decibels passband ripple
            
        elseif strcmp(FF_Type1,'Chebyshev II')
            
            [FF_N{i},FF_D{i}] = cheby2(FF_Order(i),60,1/OSR);   % 60 decibels passband ripple
        end
        
    elseif Est_FF(2,i) && Est_FF(1,i)
        
        n = FF_Order(i);
        f = [0,1/OSR,1.1/OSR,1];
        Mag = [1,1,FF_Limit(i),FF_Limit(i)];
        d = fdesign.arbmag('Nb,Na,F,A',n,n,f,Mag);
        H = design(d);
        
        H1 = H.sosMatrix;
        H2 = H.ScaleValues;
        X_N = [1];
        X_D = [1];
        H1
        for j = 1:ceil(n/2)
            X_N = conv(X_N,H1(j,1:3));
            X_D = conv(X_D,H1(j,4:6));
        end
        S1 = sum(X_N);
        S2 = sum(X_D);
        X_N = X_N*S2/S1;
        FF_N{i} = X_N
        FF_D{i} = X_D
  

    end
    try
    IF_N{i} = FF_D{i}/FF_N{i}(1)
    IF_D{i} = FF_N{i}/FF_N{i}(1)
    
    [Resp_FF{i},freq_FF{i}]=freqz(FF_N{i},FF_D{i});
    [Resp_IF{i},freq_IF{i}]=freqz(IF_N{i},IF_D{i});
    end
    
%         FF_N{1}
%         FF_D{1}
%         IF_N{1}
%         IF_D{1}
end
end