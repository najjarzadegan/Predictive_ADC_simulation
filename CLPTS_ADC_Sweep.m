function [Amp_m,Fs,BW,SNDR,ENOB,SNDR_f,ENOB_f,SNDR_t,ENOB_t,SAT_Dig] = CLPTS_ADC_Sweep(Plot_status,...
          DAC_feature, ADC_feature, SAT_ana, Amplifier_feature,...
          SH_feature, Input_feature,In_Type,H_N,H_D,LPF_N,LPF_D)

% Predictive two step ADC with shared ADC

% The system consists of a Sample and Hold (S&H), a subtractor, an ADC, a
% DAC and a digital filter working as predictor.
tic
clc
close all
%clear all
%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
% ***** Input signal **********

Amp = Input_feature{2};              % Amplitude
fin = Input_feature{3};              % frequency
N_fft = SH_feature(2);         % # of FFT points
R = SH_feature(3);             % # of repearitions (R*M cycles) where M is the number of cycles we need for N_fft samples.

% *****************************

% ***** Pipeline *************
Pipeline_Type = 1;
P_stage = 1;             % The number of pipeline stages
Count = P_stage;         % The last stage contributing to Out
LPF_status = 'OFF';       % Whether a LPF is used between the stages.
% LPF_inv_status = 'OFF';   % Whether the inverse of LPF is taken into account for digital output calculation.
%PZ_scaling = 'OFF';       % Scaling the poles and zeros of LPF between the pipeleine stages.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ***** S&H *******************

OSR = SH_feature(1);            % Over Sample Ratio

% *****************************
% ***** Amplifier *************

K = Amplifier_feature(1)*ones(P_stage,1);              % Gain

% *****************************
% ***** ADC *******************

ADC_res = ADC_feature(1)*ones(P_stage,1);       % ADC resolution
%ADC2_res = 10;
%ADC_res = [6;6];
%Max_INL = ADC_feature(2);        % Maximum INL (in LSB)
%Noise = 1;          % STD of noise added to quantizer levels

% *****************************
% ***** DAC *******************

DAC_res = DAC_feature(1)*ones(P_stage,1);       % DAC resolution
Max_INL = DAC_feature(2)*ones(P_stage,1);        % Maximum INL (in LSB)
%Noise = 1;          % STD of noise added to quantizer levels

% *****************************
% ***** Saturation *************

Sat_lim = SAT_ana*ones(P_stage,1);        % Saturation limit (dual polarity)
SAT_Dig = 2.^(ADC_res - 1);

%******************************
% ***** Predictor *************

% Everything is about the function 1/(1+kH(z)z^(-1))
%DC_gain = 1/50;        % DC gain
Order = length(H_N{1}) - 1;          % order of filter
Order_LPF = length(LPF_N{1}) - 1;  % Order of the LPF filters used between the pipeline stages
% a = 0.95*ones(P_stage,1);
% b = 0.6*ones(P_stage,1);
% theta = pi/OSR/1.0;
% Type = 1*ones(P_stage,1);
 In_type = In_Type;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%% Filter Generation %%%%%%%%%%%%%%%%%%%%%%
% L1 = zeros(1,P_stage);
% L2 = zeros(1,P_stage);
% 
% H_N = {};
% H_D = {};
% LPF_N = {};
% LPF_D = {};
% LPF_inv_N = {};
% LPF_inv_D = {};
% 
% for PL = 1:P_stage
%     G = 1;
%     
%     if Type(PL) == 1
%         %%%%%%  Type = 1  %%%%%%%%%
%         
%         p = b(PL)*ones(Order(PL),1);
%         z = a(PL)*ones(Order(PL),1);
%         
%     elseif Type(PL) == 2
%         %%%%%%  Type = 2  %%%%%%%%%
%         
%         z = 0*ones(Order(PL),1);
%         p = b(PL)*ones(Order(PL),1);
%         
%         if mod(Order(PL),2) == 1
%             z(Order(PL)) = a;
%         end
%         
%         NZ = floor(Order(PL)/2);
%         if NZ > 0
%             for i = 1:NZ
%                 
%                 A = roots([1 -2*a(PL)*cos(i*theta/NZ) a(PL)^2]);
%                 z(2*i-1) = A(1);
%                 z(2*i) = A(2);
%             end
%         end
%         
%     elseif Type(PL) == 3
%         %%%%%%  Type = 3  %%%%%%%%%
%         
%         z = 0*ones(Order(PL),1);
%         p = 0*ones(Order(PL),1);
%         
%         if mod(Order(PL),2) == 1
%             z(Order(PL)) = a(PL);
%             p(Order(PL)) = b(PL);
%         end
%         
%         NZ = floor(Order(PL)/2);
%         if NZ > 0
%             for i = 1:NZ
%                 
%                 A = roots([1 -2*a(PL)*cos(i*theta/NZ) a(PL)^2]);
%                 z(2*i-1) = A(1);
%                 z(2*i) = A(2);
%                 
%                 A = roots([1 -2*b(PL)*cos(i*theta/NZ) b(PL)^2]);
%                 p(2*i-1) = A(1);
%                 p(2*i) = A(2);
%                 
%             end
%         end
%         
%         
%     end
%     
%     [N,D] = zp2tf(z,p,G);
%     
%     %%%%%%%%%%%%%%
%     
%     % N = [0.079 -0.454 1.079 -1.366 0.973 -0.37 0.058];
%     % D = [1 -0.993 0.905 -0.277 0.114 0.00158 0.0049];
%     % Order = length(N) - 1
%     
%     %%%%%%%%%%%%%
%     
%     H_N1 = (D - N)/K(PL);
%     H_D1 = N;
%     
%     H_N{PL} = H_N1;
%     H_D{PL} = H_D1;
%     
%     N_noise = K(PL)*conv(N,H_N1);
%     D_noise = conv(D,H_D1);
%     
%     try
%         L2(PL) = filternorm(N_noise,D_noise);
%     catch ME
%         warning('The filter is almost unstable')
%     end
%     
%     [h,w]=freqz(N,D);
%     L1(PL) = max(abs(h(1:length(w)/OSR)));
%     
%     
%     [LPF_N{PL},LPF_D{PL}] = butter(Order_LPF(PL),1/OSR,'low');
%     
%     if strcmp(PZ_scaling,'ON') && PL < P_stage
%         [z,p,k] = butter(Order_LPF(PL),1/OSR,'low');
%         Max = max(max(abs(z),abs(p)));
%         CF = 0.95/Max;
%         
%         LPF_N{PL} = LPF_N{PL}.*CF.^(0:Order_LPF(PL));
%         %roots(LPF_N{PL})
%         LPF_D{PL} = LPF_D{PL}.*CF.^(0:Order_LPF(PL));
%         %roots(LPF_D{PL})
%     end
%     
%     LPF_inv_D{PL} = LPF_N{PL}/LPF_N{PL}(1);
%     LPF_inv_N{PL} = LPF_D{PL}/LPF_N{PL}(1);
%     
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENOB(1) = 0.6*(ADC_res(P_stage) + sum(log2(K)));
j = 2;
A(1) = 0.0;

figure
xlabel('Amplitude')
ylabel('ENOB')
title('ENOB vs Amplitude')
h = animatedline;
axis([0 1.2 0 (ADC_res(P_stage) + sum(log2(K))) + 1])
addpoints(h,A(1),ENOB(1))
drawnow

while (ENOB(j-1) > 0.5*(ADC_res(P_stage) + sum(log2(K))) || j < 10) 
    
    
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
Fs = 1/Ts_OS;
BW = Fs/2/OSR;

step = Ts_OS/10;       % Time step
t = step:step:P;       % Time vector (number of points is N_fft*R*OSR*10)

if strcmp(In_type,'Sine')
    X = A(j)*sin(2*pi*fin*t); % Input Signal
elseif strcmp(In_type,'Chirp')
    X = A(j)*chirp(t,0,P,fin); % Input Signal
elseif strcmp(In_type,'Multitone')
    X = sin(2*pi*fin*t) + sin(2*pi*fin/2*t) + sin(2*pi*fin/4*t) + sin(2*pi*3*fin/4*t); % Input Signal
    X = A(j)*X/max(X);
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
DAC_i = zeros(P_stage,L);         % Initializing DAC
Out = zeros(1,L);               % Output digital
Out_f = zeros(1,L);             % Filtered Out
Out_bf = zeros(P_stage,L);  % A vector storing the digital number Out so that applies it to the inverse filter
Err = zeros(1,L);               % Whole ADC Error
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
    %%%%% Pipelining: LPF between the stages %%%%%%%%%%
    
    if  i < L
        for p = 2:P_stage
            
            if strcmp(LPF_status,'ON')
                
                %%%%%%%%%%%%%%%%%%%  LPF between the stages    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if i < Order_LPF(p - 1) + 1
                    S2(p,i + 1) = 0;
                else
                    
                    S2(p,i + 1) = S2(p,i - Order_LPF(p - 1) + 1:i)*(flip(-LPF_D{p - 1}(2:end)))' + ...
                        E2(p - 1,i - Order_LPF(p - 1):i)*(flip(LPF_N{p - 1}(1:end)))';
                    
                    %Pred(i) = Pred(i) + H_N(1)*(2*ADC(i) - ADC(i -1));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            else
                S2(p,i + 1) = E2(p - 1,i);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
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
%             Pred(p,i - Order(p):i - 1)
%             ADC(p,i - Order(p) + 1:i)
            
            Pred(p,i) = Pred(p,i - Order(p):i - 1)*(flip(-H_D{p}(2:end)))' + ...
                ADC(p,i - Order(p) + 1:i)*(flip(H_N{p}(2:end)))';
            
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
    Y3 = floor(Y1) + 0.5;
    
    if i < L
        DAC(:,i+1)   = Sat_lim.*Y2./2.^(DAC_res - 1);
        DAC_i(:,i+1) = Sat_lim.*Y3./2.^(DAC_res - 1);
        
    end
    
    
end

%%%%%%%%%%%%%%%%%% Pipeline Type 2 %%%%%%%%%%%%%%%%%%%%%%%%


% In this topology, the first stage is the proposed ADC (Closed-loop
% Two-step) and the second stage is a conventinal ADC.

if Pipeline_Type == 2

ADC2 = floor(2^(ADC2_res - 1)*SAT(1,:)/Sat_lim);
%ADC2 = E2(1,:);

Out = (ADC2 + 0.5)/(2^(ADC2_res - 1))/K(1) + DAC_i(1,:);
%Out = (ADC2)/K(1) + DAC(1,:);


%%%%%%%%%%%%%%%%%%  Digital domain processing %%%%%%%%%%%%%%%%%

%     if i - (Count -1) <1
%         Err(i) = 0;
%         Out(i) = 0;
%     else

elseif Pipeline_Type == 1

Out = (ADC(Count,:) + 0.5)/(2^(ADC_res(Count) - 1));

for pp = 1:Count
    p = Count - pp + 1;
    Out = Out/K(p) + DAC_i(p,:);
    %(p - Count)
    % + (p - Count)
    Out_bf(p,:) = Out;
    
    %%%%%%%%%%%%%%% Inverse Filter %%%%%%%%%%%%%%%%%%
%     if strcmp(LPF_inv_status,'ON')
%         
%         if p > 1
%             
%             for i = 1:L-1
%                 
%                 if i < Order_LPF(p - 1) + 1
%                     Out(i) = E2(p - 1,i);
%                 else
%                     
%                     Out(i) = Out(i - Order_LPF(p - 1):i - 1)*(flip(-LPF_inv_D{p - 1}(2:end)))' + ...
%                         Out_bf(p,i - Order_LPF(p - 1) + 1:i + 1)*(flip(LPF_inv_N{p - 1}(1:end)))';
%                     
%                 end
%                 
%                 %                     S_O(i+1) = S2(2,i) - Out_bf(2,i);
%                 %                     if i < Order_LPF(p - 1) + 1
%                 %                         E2_est(i) = 0;
%                 %                     else
%                 %
%                 %                         E2_est(i) = E2_est(i - Order_LPF(p - 1):i - 1)*(flip(-LPF_inv_D{p - 1}(2:end)))' + ...
%                 %                             S_O(i - Order_LPF(p - 1)+1:i+1)*(flip(LPF_inv_N{p - 1}(1:end)))';
%                 %
%                 %                     end
%             end
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
end
%Out(i) = (ADC(2,i) + 0.5)/(2^(ADC_res(2) - 1))/K(1)/K(2) + DAC(2,i)/K(1) + DAC(1,i-1);

%         if i - (Count -1) <1
%             Err(i) = 0;
%             %Out(i) = 0;
%         else
%             Err(i) = Out(i) - S2(1,i-Count+1);
%         end
Err = Out - S2(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Out2 = zeros(1,L);
%   Err2 = zeros(1,L);
%             Out2(i) = (ADC(2,i) + 0.5)/(2^(ADC_res(2) - 1))/K(2) + DAC(2,i);
%
%             Err2(i) = Out2(i) - S2(2,i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%  Output Filter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:L - 1
    if i < Order_LPF(P_stage) + 1
        Out_f(i + 1) = 0;
    else
        
        Out_f(i + 1) = Out_f(i - Order_LPF(P_stage) + 1:i)*(flip(-LPF_D{P_stage}(2:end)))' + ...
            Out(i - Order_LPF(p):i)*(flip(LPF_N{P_stage}(1:end)))';
        
    end
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%   NEW METHOD FOR SNR CALCULATION %%%%%%%%%%%%%%%%%%%%%%
% 
% % % Max_out = max(Out(10*N_fft*OSR:14*N_fft*OSR));
% % % Ph_out = find(Out(10*N_fft*OSR:14*N_fft*OSR) == Max_out);
% %
% Max_out_f = max(Out_f(10*N_fft*OSR:14*N_fft*OSR));
% Ph_out_f = find(Out_f(10*N_fft*OSR:14*N_fft*OSR) == Max_out_f);
% 
% Max_XS = max(X_S(10*N_fft*OSR:14*N_fft*OSR));
% Ph_XS = find(X_S(10*N_fft*OSR:14*N_fft*OSR) == Max_XS);
% 
% %r1 = Ph_out(1) - Ph_XS(1);
% XS_new = X_S(10*N_fft*OSR:50*N_fft*OSR);
% 
% %Out_new = max(X_S)*Out(10*N_fft*OSR + r1:50*N_fft*OSR + r1)/Max_out;
% 
% r2 = Ph_out_f(1) - Ph_XS(1);
% 
% for i =1:11
%     sh = i - 6;
%     for j = 1:11
%         G = 1 - (j-6)/50;
%         for k = 1:11
%             DC = (k-6)/100;
%             Out_f_new = G*max(X_S)*Out_f(10*N_fft*OSR + r2 + sh:50*N_fft*OSR + r2 + sh)/Max_out_f + DC;
%             
%             % plot(10*N_fft*OSR:50*N_fft*OSR,Out_new,10*N_fft*OSR:50*N_fft*OSR,XS_new)
%             % legend('Out_new','X_S')
%             %
%             %             figure
%             %             plot(10*N_fft*OSR:50*N_fft*OSR,Out_f_new,10*N_fft*OSR:50*N_fft*OSR,XS_new)
%             %             legend('Out_f_new','X_S')
%             
%             Err_new = Out_f_new - XS_new;
%             SNR_new(i,j,k) = 20*log10(norm(XS_new)/norm(Err_new));
%             ENOB_new(i,j,k) = (SNR_new(i,j,k) - 1.76)/6.02;
%             [ENOB_new(i,j,k),Amp];
%             
%         end
%     end
% end
% MM = max(max(max(ENOB_new)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sat_status = zeros(1,P_stage);
% for p = 1:P_stage
%     Sat_status(p) = max(abs(E2(p,10*N_fft*OSR:L)));
% end
%norm(Err(10*N_fft*OSR:L))

SNDR_t = 20*log10(norm(S2(1,10*N_fft*OSR:L))/norm(Err(10*N_fft*OSR:L)));
%     SNR2 = 20*log10(norm(S2(2,10*N_fft*OSR:L))/norm(Err2(10*N_fft*OSR:L)));
%norm(S2(N_fft*OSR:L))
%norm(Err(N_fft*OSR:L))

ENOB(j) = (SNDR_t - 1.76)/6.02;
%     ENOB2(j) = (SNR2 - 1.76)/6.02;
% [ENOB,Amp]
%     [ENOB2(j),A(j)]
%ENOB2 = ENOB + log2(Sat_lim/A);

if ENOB(j) > max(ENOB(1:j - 1))

Out_m = Out;
X_S_m = X_S;
Out_f_m = Out_f;
E2_m = E2;
Pred_m = Pred;
X_m = X;
S1_m = S1;
SAT2_m = SAT2;
SAT_m =SAT;
DAC_m = DAC;
ADC_m = ADC;
Amp_m = A(j);
ENOB_t = ENOB(j);

end

addpoints(h,A(j),ENOB(j))
drawnow

j = j + 1;

end


% legend(['Gain = ' num2str(K),', OSR = ' num2str(OSR),...
%     ADC Res = ' num2str(ADC_res),', DAC Res = ' num2str(DAC_res),...
%     Type = ' num2str(Type),', Order = ' num2str(Order),...
%     Alpha = ' num2str(a),', Beta = ' num2str(b)]);

% Term2 = L1
% Term3 = L2
% M0 = log2((1 + L2)./(1./K' - max(Amp)*L1))

% plot(t,X,t,S1,'-')
% legend('Input', 'Sampled and Held')
%
% figure
% plot(E2)
% legend('Amplified')
%
if strcmp(Plot_status{9},'Only Out')
    
    try
        figure
        plot(t(10:10:end),Out_m)
        title('Time domain output')
        xlabel('Time')
        %legend('Output Digital')
    end
elseif strcmp(Plot_status{9},'Out and Input')
    try
        figure
        plot(t(10:10:end),Out_m,t(10:10:end),X_S_m)
        title('Time domain input and output')
        xlabel('Time')
        legend('Output Digital','Input Signal')
    end
end
%
if strcmp(Plot_status{11},'Only Out')
    try
        figure
        plot(t(10:10:end),Out_f_m)
        title('Time domain filtered output')
        xlabel('Time')
        %legend('Filtered Output')
    end
elseif strcmp(Plot_status{11},'Out and Input')
    try
        figure
        plot(t(10:10:end),Out_f_m,t(10:10:end),X_S_m)
        title('Time domain filtered output and Input')
        xlabel('Time')
        legend('Filtered Output','Input Signal')
    end
end
%
if Plot_status{10}
    try
        figure
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(Out_m(L/2 + 1:L)))/max(abs(fft(Out_m(L/2 + 1:L))))))
        ylabel('dB')
        xlabel('Frequency')
        title('Digital Out FFT')
        %     L_Out_f = log2(length(Out_f))
    end
end
%
if Plot_status{12}
    try
        figure
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(Out_f_m(L/2 + 1:L)))/max(abs(fft(Out_f_m(L/2 + 1:L))))))
        ylabel('dB')
        xlabel('Frequency')
        title('Digital Out Filtered FFT')
        %     L_Out_f = log2(length(Out_f))
    end
end
%
if Plot_status{13}
    try
        figure
        semilogy(abs(fft(E2_m(1,L/2 + 1:L))))
        title('Amplifier output frequency domain')
        xlabel('Frequency')
        legend('Amplifier FFT')
        %     L_E2 = log2(length(E2(1,:)))
    end
end
%
if Plot_status{4}
    try
        figure
        plot(t(10:10:end),E2_m(1,:))
        title('Amplifier output time domain')
        xlabel('Time')
        legend('Amplifier output')
    end
end
%
if Plot_status{1}
    try
        figure
        plot(t(10:10:end),Pred_m(1,:))
        xlabel('Time')
        title('Prediction')
    end
end
%
if Plot_status{2} && Plot_status{3}
    try
        figure
        plot(t,X_m,t,S1_m)
        xlabel('Time')
        title('Input ans Sampled Signals')
        legend('Input', 'Sampled')
    end
elseif Plot_status{2}
    try
        figure
        plot(t,X_m)
        xlabel('Time')
        title('Input Signals')
    end
elseif Plot_status{3}
    try
        figure
        plot(t,S1_m)
        xlabel('Time')
        title('Sampled Signals')
    end
end
%
if Plot_status{5}
    try
        figure
        plot(t(10:10:end),SAT_m(1,:))
        xlabel('Time')
        title('Analog Saturation block')
    end
end
%
if Plot_status{6}
    try
        figure
        plot(t(10:10:end),SAT2_m(1,:))
        xlabel('Time')
        title('Digital Saturation block')
    end
end
%
if Plot_status{7}
    try
        figure
        plot(t(10:10:end),ADC_m(1,:))
        xlabel('Time')
        ylabel('Digital Code')
        title('ADC output')
    end
end
%
if Plot_status{8}
    try
        figure
        plot(t(10:10:end),DAC_m(1,:))
        xlabel('Time')
        title('DAC output')
    end
end
%
if Plot_status{14}
    try
        figure
        plot(t(10:10:end),X_S_m,t(10:10:end),DAC_m)
        xlabel('Time')
        title('Prediction vs actual input')
        legend('Actual input','Prediction')
    end
end
% try
%     figure
%     semilogy(abs(fft(S2(2,:))))
%     legend('S2 frequency response')
% %     L_S2 = log2(length(S2(1,:)))
% end
%
% figure
% plot(ADC)
% legend('ADC')
%
% figure
% plot(DAC)
% hold on
% plot(S2)
% legend('DAC','Sample')
%
%
% figure
% semilogy(abs(fft(Out)))
% legend('Digital Out')
% %
% try
%     figure
%     plot(20*log10(abs(fft(Out))/max(abs(fft(Out)))))
%     ylabel('dB')
%     legend('Digital Out')
% %     L_Out = log2(length(Out))
% end
% %
% figure
% semilogy(abs(fft(Out_f)))
% legend('Digital Out Filtered')
% %
% try
%     figure
%     plot(20*log10(abs(fft(Out_f))/max(abs(fft(Out_f)))))
%     ylabel('dB')
%     legend('Digital Out Filtered')
% %     L_Out_f = log2(length(Out_f))
% end
%
% try
%     figure
%     plot(1:P_stage,Sat_status)
% end
%
% try
%     figure
%     plot(1:length(Out_bf(2,:)),Out_bf(2,:),1:length(Out_bf(2,:)),S2(2,:))
%     legend('Out_bf','S2(2)')
%     max(abs(Out_bf(2,100:end)-S2(2,100:end)))
% end
%
% figure
% plot(1:length(E2_est),E2_est,1:length(E2_est),E2(1,:))
% legend('E2_est','E2(1)')
% max(abs(E2_est - E2(1,:)))

% try
%     figure
%     plot(1:length(Out),Out,1:length(Out),E2(1,:))
%     legend('Out','E2(1)')
%     max(abs(Out - E2(1,:)))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFT_D_O = abs(fft(Out_m(L/2 + 1:L)));
%  figure
%  plot(FFT_D_O)
Max_bin = max(FFT_D_O(1:L/4));
Fm = find(FFT_D_O(1:L/4) == Max_bin);

FFT_D_O_new = FFT_D_O(1:L/4);
FFT_D_O_new(Fm) = 0;

% Max_bin^2
% sum(FFT_D_O_new.^2)

SNDR = 10*log10(Max_bin^2/sum(FFT_D_O_new.^2));
ENOB = (SNDR - 1.76)/6.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FFT_D_O = abs(fft(Out_f_m(L/2 + 1:L)));
%  figure
%  plot(FFT_D_O)
Max_bin = max(FFT_D_O(1:L/4));
Fm = find(FFT_D_O(1:L/4) == Max_bin);

FFT_D_O_new = FFT_D_O(1:L/4);
FFT_D_O_new(Fm) = 0;

%Max_bin^2
%sum(FFT_D_O_new.^2)

SNDR_f = 10*log10(Max_bin^2/sum(FFT_D_O_new.^2));
ENOB_f = (SNDR_f - 1.76)/6.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
