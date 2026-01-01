% Predictive two step ADC with shared ADC

% The system consists of a Sample and Hold (S&H), a subtractor, an ADC, a
% DAC and a digital filter working as predictor.
tic
close all
clear all
%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
% ***** Input signal **********

%A = 0.5;              % Amplitude
fin = 1;            % frequency
N_fft = 128;         % # of FFT points
R = 64;             % # of repearitions (R*M cycles) where M is the number of cycles we need for N_fft samples.

% *****************************

% ***** Pipeline *************
P_stage = 1;        % The number of pipeline stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ***** S&H *******************

OSR = 16;            % Over Sample Ratio

% *****************************
% ***** Amplifier *************

K = 16*ones(P_stage,1);              % Gain

% *****************************
% ***** Saturation *************

Sat_lim = 1*ones(P_stage,1);        % Saturation limit (dual polarity)

%******************************
% ***** ADC *******************

ADC_res = 6*ones(P_stage,1);       % ADC resolution
%Max_INL = 1;        % Maximum INL (in LSB)
%Noise = 1;          % STD of noise added to quantizer levels

% *****************************
% ***** DAC *******************

DAC_res = 6*ones(P_stage,1);       % DAC resolution
Max_INL = 0*ones(P_stage,1);        % Maximum INL (in LSB)
%Noise = 1;          % STD of noise added to quantizer levels

% *****************************
% ***** Predictor *************

% Everything is about the function 1/(1+kH(z)z^(-1))
%DC_gain = 1/50;        % DC gain
Order = 2*ones(P_stage,1);          % order of filter
a = 0.95*ones(P_stage,1);
b = 0.0*ones(P_stage,1);
theta = pi/OSR/1.0;
Type = 1*ones(P_stage,1);
In_type = 'sin';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Filter Generation %%%%%%%%%%%%%%%%%%%%%%
L1 = zeros(1,P_stage);
L2 = zeros(1,P_stage);

H_N = {};
H_D = {};

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
    
    H_N1 = (D - N)/K(PL)
    H_D1 = N
    
    H_N{PL} = H_N1;
    H_D{PL} = H_D1;
    
    N_noise = K(PL)*conv(N,H_N1);
    D_noise = conv(D,H_D1);
    
    try
        L2(PL) = filternorm(N_noise,D_noise);
    catch ME
        warning('The filter is almost unstable')
    end
    
    [h,w]=freqz(N,D);
    L1(PL) = max(abs(h(1:length(w)/OSR)));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENOB(1) = 0.6*(ADC_res(P_stage) + sum(log2(K)));
j = 2;
A(1) = 0.0;
while (ENOB(j-1) > 0.5*(ADC_res(P_stage) + sum(log2(K))) || j < 15) 
    
    
    if ENOB(j - 1) < (ADC_res(P_stage) + sum(log2(K))) - 2
        A(j) = 0.1 + A(j - 1);
    elseif (ENOB(j - 1) > (ADC_res(P_stage) + sum(log2(K))) - 2) && (ENOB(j - 1) < (ADC_res(P_stage) + sum(log2(K))) - 1)
        A(j) = 0.05 + A(j - 1);
    else
        A(j) = 0.02 + A(j - 1);
    end
    %%%%%%%%%%%% Parts %%%%%%%%%%%%%%%%%%%%
    
    % ***** Input Signal *******************
    
    M = N_fft/2 - 1;       % The number of cycles we need for N_fft samples
    Tin = 1/fin;           % Input period;
    P = Tin*M*R;           % Simulation period
    
    Ts_Nyq = Tin*M/N_fft;  % Sampling Period for nyquist rate sampling
    Ts_OS = Ts_Nyq/OSR ;   % Sampling Period for oversampling
    
    step = Ts_OS/10;       % Time step
    t = step:step:P;       % Time vector (number of points is N_fft*R*OSR*10)
    
    if strcmp(In_type,'sin')
        X = A(j)*sin(2*pi*fin*t); % Input Signal
    elseif strcmp(In_type,'chirp')
        X = A(j)*chirp(t,0,P,fin/1); % Input Signal
    end
    
    % ***** Sample and Hold ***************
    
    S1 = X;
    for i = 1:N_fft*R*OSR*10
        
        if (mod(i,10) > 5)
            S1(i) = X(10*round(i/10) - 5);
        end
    end
    
    % Sampled and Held input
    
    X_S = X(10:10:N_fft*R*OSR*10);
    L = length(X_S);
    S2 = zeros(P_stage,L);
    S2(1,:) = X_S;
    % Sampled input
    
    
    
    %%%%%% Feedback Loop %%%%%%%%%%%%%%%%%%
    
    % # of samples
    E1 = zeros(P_stage,L);          % Loop Error
    E2 = zeros(P_stage,L);          % Error after amplification
    SAT = zeros(P_stage,L);         % Saturation block output
    ADC = zeros(P_stage,L);         % ADC output
    Pred = zeros(P_stage,L);        % Predictor output
    SAT2 = zeros(P_stage,L);        % saturation output
    DAC = zeros(P_stage,L);         % Initializing DAC
    Out = zeros(1,L);         % Output digital
    Err = zeros(1,L);         % Whole ADC Error
    G = zeros(P_stage,L);           % Adaptive gain
    
    for i = 1:L
        
        %%%%% Subtractor  %%%%%%%%%
        %i;
        E1(:,i) = S2(:,i) - DAC(:,i);
        %S2(i);
        %S2(i) - DAC(i);
        %%%%% Gain stage %%%%%%%%%%
        
        %G(i) = (1 - exp(-0.0001*i))*K;
        E2(:,i) = K.*E1(:,i);
        %K*E1(i)
        %%%%% Pipelining %%%%%%%%%%
        
        for p = 2:P_stage
            S2(p,i + 1) = E2(p - 1,i);
        end
        
        %%%%% Saturation %%%%%%%%%%
        
        SAT(:,i) = E2(:,i);
        
        for p = 1:P_stage
            if (E2(p,i) > Sat_lim(p))
                SAT(p,i) = Sat_lim(p);
                
            elseif (E2(p,i) < -Sat_lim(p))
                SAT(p,i) = -Sat_lim(p);
                
            end
        end
        
        %%%% ADC %%%%%%%%%%%%%%%%%%
        
        ADC(:,i) = floor(2.^(ADC_res - 1).*SAT(:,i)./Sat_lim);
        %floor(2^(ADC_res - 1)*SAT(i)/Sat_lim)
        % Nonidealities like INL and noise should be added later
        
        %%%% Predictor %%%%%%%%%%%%
        %
        % Type 1: zeros close to unit circle and poles in center
        % (Z^order - (1 - DC_gain))/Z^order
        %     if Type == 1
        %         if (i - Order < 1)
        %             %Pred(i) = 2^(ADC_res - 1)*S2(i);
        %             Pred(i) = 0;
        %          %   0
        %         else       Pred(i) = (a)*Pred(i - Order) + ADC(i - Order + 1)*(a - b)/K;
        %           %  (1 - DC_gain)*Pred(i - Order) + ADC(i - Order)*(1 - DC_gain)/K
        %         end
        % %     Pred(1) = ADC(1);
        % %     Pred(2) = ADC(2);
        % %
        % %     if (i>2)
        % %         Pred(i) = 2*(1 - DC_gain)*Pred(i - 2) - (1 - DC_gain)^2*Pred(i - 1)...
        % %         - ADC(i - 1)*2*(1 - DC_gain)/K + (1 - DC_gain)^2*Pred(i - 2)/K;
        %      end
        
        %
        %     % Type 2: zeros close to 1 and poles in center
        %     %
        %     elseif Type == 2
        %         if (i - 2 < 1)
        %             %Pred(i) = 2^(ADC_res - 1)*S2(i);
        %             Pred(i) = 0;
        %            %0
        %         else
        %             Pred(i) = 2*a*cos(theta)*Pred(i - 1) - a^2*Pred(i - 2) + (-a^2/K)*ADC(i - 1) + 2*a*cos(theta)*ADC(i)/K;
        %            %(1 - DC_gain)*Pred(i - Order) + ADC(i - Order)*(1 - DC_gain)/K
        %         end
        %
        %     elseif Type == 3
        %
        %         % Type 3: 3 zeros close to 1 and poles located at 'b'.
        %         if (i < 4)
        %             %Pred(i) = 2^(ADC_res - 1)*S2(i);
        %             Pred(i) = 0;
        %            %0
        %         else
        %             Pred(i) = (a + 2*a*cos(theta))*Pred(i - 1) + (-a^2 - 2*a^2*cos(theta))*Pred(i - 2) + a^3*Pred(i - 3) + ...
        %             (a^3 - b^3)*ADC(i - 2)/K + (3*b^2 - a^2 - 2*a^2*cos(theta))*ADC(i - 1)/K + (-3*b + a + 2*a*cos(theta))*ADC(i)/K;
        %            %(1 - DC_gain)*Pred(i - Order) + ADC(i - Order)*(1 - DC_gain)/K
        %         end
        %     end
        
        %%%%%%%%%%%%%%%%%%%  NEW Predictor   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for p = 1:P_stage
            
            if i < Order(p) + 1
                Pred(p,i) = 0;
            else
                %         Out(i - Order:i - 1)
                %         (flip(-H_D(2:end)))'
                %         In(i - Order + 1:i)
                %         (flip(H_N(2:end)))'
                Pred(p,i) = Pred(p,[i - Order(p):i - 1])*(flip(-H_D{p}(2:end)))' + ...
                    ADC(p,[i - Order(p) + 1:i])*(flip(H_N{p}(2:end)))';
                
                %Pred(i) = Pred(i) + H_N(1)*(2*ADC(i) - ADC(i -1));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%% Saturation block %%%%%
        
        SAT2(:,i) = Pred(:,i);
        
        for p = 1:P_stage
            
            if (Pred(p,i) > 2^(ADC_res(p) - 1))
                SAT2(p,i) = 2^(ADC_res(p) - 1);
            elseif (Pred(p,i) < -2^(ADC_res(p) - 1))
                SAT2(p,i) = -2^(ADC_res(p) - 1);
            end
        end
        
        %%%% DAC %%%%%%%%%%%%%%%%%%
        
        Y1 = SAT2(:,i)./2.^(ADC_res - DAC_res);
        Y2 = floor(Y1) + 0.5 - Max_INL.*sin(pi*Y1./2.^(DAC_res - 1));
        
        DAC(:,i+1) = Sat_lim.*Y2./2.^(DAC_res - 1);
        %Sat_lim*Y2/2^(DAC_res - 1)
        
        Count = P_stage;        % The last stage contributing to Out
        
        if i - (Count -1) <1
            Err(i) = 0;
            Out(i) = 0;
        else
            
            Out(i) = (ADC(Count,i) + 0.5)/(2^(ADC_res(Count) - 1));
            for pp = 1:Count
                p = Count - pp + 1;
                Out(i) = Out(i)/K(p) + DAC(p,i + (p - Count));
            end
            %Out(i) = (ADC(2,i) + 0.5)/(2^(ADC_res(2) - 1))/K(1)/K(2) + DAC(2,i)/K(1) + DAC(1,i-1);
            
            
            Err(i) = Out(i) - S2(1,i - (Count -1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Out2 = zeros(1,L);
        %   Err2 = zeros(1,L);
        %             Out2(i) = (ADC(2,i) + 0.5)/(2^(ADC_res(2) - 1))/K(2) + DAC(2,i);
        %
        %             Err2(i) = Out2(i) - S2(2,i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    end
    
    for p = 1:P_stage
        Sat_status(p,j) = max(abs(E2(p,10*N_fft*OSR:L)));
    end
    %norm(Err(10*N_fft*OSR:L))
    
    SNR = 20*log10(norm(S2(1,10*N_fft*OSR:L))/norm(Err(10*N_fft*OSR:L)));
    %     SNR2 = 20*log10(norm(S2(2,10*N_fft*OSR:L))/norm(Err2(10*N_fft*OSR:L)));
    %norm(S2(N_fft*OSR:L))
    %norm(Err(N_fft*OSR:L))
    
    ENOB(j) = (SNR - 1.76)/6.02;
    %     ENOB2(j) = (SNR2 - 1.76)/6.02;
    [ENOB(j),A(j)]
    %     [ENOB2(j),A(j)]
    %ENOB2 = ENOB + log2(Sat_lim/A);
    j = j +1;
end

Peak_ENOB = max(ENOB)
Peak_SNR = 6.02*Peak_ENOB + 1.76
Amp = A(find(ENOB==Peak_ENOB))


plot(A,Sat_status)

figure
plot(A,ENOB)
ylabel('ENOB')
xlabel('Amplitude')
% legend(['Gain = ' num2str(K),', OSR = ' num2str(OSR),...
%     ADC Res = ' num2str(ADC_res),', DAC Res = ' num2str(DAC_res),...
%     Type = ' num2str(Type),', Order = ' num2str(Order),...
%     Alpha = ' num2str(a),', Beta = ' num2str(b)]);

Term2 = L1
Term3 = L2
M0 = log2((1 + L2)./(1./K - max(A)*L1))

% plot(t,X,t,S1,'-')
% legend('Input', 'Sampled and Held')
%
% figure
% plot(E2)
% legend('Amplified')
%
% figure
% plot(ADC)
% legend('ADC')
%
% figure
% plot(DAC)
%
% hold on
%
% plot(S2)
% legend('DAC','Sample')

toc
