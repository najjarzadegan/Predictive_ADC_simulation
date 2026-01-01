function [theta,H_N,H_D,LPF_N,LPF_D,h_loop,w_loop,h_H,w_H,h_LPF,w_LPF,Term] = ...
    PredictorFilter(K,Predictor_feature,P_Type,LPF_order,LPF_Type)

P_stage = 1;
PZ_scaling = 'OFF';

OSR = Predictor_feature(5);
Order = Predictor_feature(3)*ones(P_stage,1);          % order of filter
Order_LPF = LPF_order*ones(P_stage,1);  % Order of the LPF filters used between the pipeline stages
a = Predictor_feature(1)*ones(P_stage,1);
b = Predictor_feature(2)*ones(P_stage,1);
r = Predictor_feature(4);
theta = 0;
if P_Type > 1
    theta = pi/OSR/r;
end
Type = P_Type*ones(P_stage,1);
%In_type = 'sin';

%%%%%%%% Filter Generation %%%%%%%%%%%%%%%%%%%%%%
L1 = zeros(1,P_stage);
L2 = zeros(1,P_stage);


H_N = {};
H_D = {};
LPF_N = {};
LPF_D = {};
LPF_inv_N = {};
LPF_inv_D = {};

for PL = 1:P_stage
    G = 1;
    
    if Type(PL) == 1
        %%%%%%  Type = 1  %%%%%%%%%
        
        p = b(PL)*ones(Order(PL),1);
        z = a(PL)*ones(Order(PL),1);
        
    elseif Type(PL) == 2
        %%%%%%  Type = 2  %%%%%%%%%
        
        z = 0*ones(Order(PL),1);
        p = b(PL)*ones(Order(PL),1);
        
        if mod(Order(PL),2) == 1
            z(Order(PL)) = a;
        end
        
        NZ = floor(Order(PL)/2);
        if NZ > 0
            for i = 1:NZ
                
                A = roots([1 -2*a(PL)*cos(i*theta/NZ) a(PL)^2]);
                z(2*i-1) = A(1);
                z(2*i) = A(2);
            end
        end
        
    elseif Type(PL) == 3
        %%%%%%  Type = 3  %%%%%%%%%
        
        z = 0*ones(Order(PL),1);
        p = 0*ones(Order(PL),1);
        
        if mod(Order(PL),2) == 1
            z(Order(PL)) = a(PL);
            p(Order(PL)) = b(PL);
        end
        
        NZ = floor(Order(PL)/2);
        if NZ > 0
            for i = 1:NZ
                
                A = roots([1 -2*a(PL)*cos(i*theta/NZ) a(PL)^2]);
                z(2*i-1) = A(1);
                z(2*i) = A(2);
                
                A = roots([1 -2*b(PL)*cos(i*theta/NZ) b(PL)^2]);
                p(2*i-1) = A(1);
                p(2*i) = A(2);
                
            end
        end
        
        
    end
    
    [N,D] = zp2tf(z,p,G);
    
    %%%%%%%%%%%%%%
    
    % N = [0.079 -0.454 1.079 -1.366 0.973 -0.37 0.058];
    % D = [1 -0.993 0.905 -0.277 0.114 0.00158 0.0049];
    % Order = length(N) - 1
    
    %%%%%%%%%%%%%
    
    H_N1 = (D - N)/K(PL);
    H_D1 = N;
    
    H_N{PL} = H_N1;
    H_D{PL} = H_D1;
    
    N_noise = K(PL)*conv(N,H_N1);
    D_noise = conv(D,H_D1);
    
    try
        L2(PL) = filternorm(N_noise,D_noise);
    catch ME
        warning('The filter is almost unstable')
    end
    
    N
    D
    [h_loop,w_loop] = freqz(N,D);
    [h_H,w_H] = freqz(H_N{PL},H_D{PL});
    
    
    L1(PL) = max(abs(h_loop(1:length(w_loop)/OSR)));
    
    if strcmp(LPF_Type,'Butterworth')
        
        [LPF_N{PL},LPF_D{PL}] = butter(Order_LPF(PL),1/OSR,'low');
        
    elseif strcmp(LPF_Type,'FIR')
        
        %[LPF_N{PL},LPF_D{PL}] = butter(Order_LPF(PL),1/OSR,'low');
        
    elseif strcmp(LPF_Type,'Chebyshev I')
        
        %[LPF_N{PL},LPF_D{PL}] = butter(Order_LPF(PL),1/OSR,'low');
        
    elseif strcmp(LPF_Type,'Chebyshev II')
        
        %[LPF_N{PL},LPF_D{PL}] = butter(Order_LPF(PL),1/OSR,'low');
        
    end
    
    [h_LPF,w_LPF]=freqz(LPF_N{PL},LPF_D{PL});
    
    if strcmp(PZ_scaling,'ON') && PL < P_stage
        [z,p,k1] = butter(Order_LPF(PL),1/OSR,'low');
        Max = max(max(abs(z),abs(p)));
        CF = 0.95/Max;
        
        LPF_N{PL} = LPF_N{PL}.*CF.^(0:Order_LPF(PL));
        %roots(LPF_N{PL})
        LPF_D{PL} = LPF_D{PL}.*CF.^(0:Order_LPF(PL));
        %roots(LPF_D{PL})
    end
    
    LPF_inv_D{PL} = LPF_N{PL}/LPF_N{PL}(1);
    LPF_inv_N{PL} = LPF_D{PL}/LPF_N{PL}(1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Term2 = L1;
Term3 = L2;
M0 = log2((1 + L2)./(1./K' - 1*L1));
DoC = (1./K')./(1./K' - 1*L1); % Degree of Confidence
Term = [L1,L2,M0,DoC];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end