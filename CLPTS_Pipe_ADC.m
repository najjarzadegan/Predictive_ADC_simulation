function [Fs,BW,SNDR,ENOB,SNDR_f,ENOB_f,SNDR_t,ENOB_t,SAT_Dig] = CLPTS_ADC(Plot_status,...
          DAC_feature, ADC_feature, SAT_ana, Amplifier_feature,...
          SH_feature, Input_feature,In_Type,H_N,H_D,LPF_N,LPF_D,Pipeline_Type,P_stage,FF_N,FF_D,IF_N,IF_D)

% Predictive two step ADC with shared ADC

% The system consists of a Sample and Hold (S&H), a subtractor, an ADC, a
% DAC and a digital filter working as predictor.
tic
clc
close all
%clear all
% FF_N{1}
% FF_D{1}
% IF_N{1}
% IF_D{1}
%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
% ***** Input signal **********
Plot_status;
Amp = Input_feature{2};              % Amplitude
fin = Input_feature{3};              % frequency
N_fft = SH_feature(2);         % # of FFT points
R = SH_feature(3);             % # of repearitions (R*M cycles) where M is the number of cycles we need for N_fft samples.

% *****************************

% ***** Pipeline *************
Pipeline_Type = Pipeline_Type;
P_stage = str2double(P_stage);             % The number of pipeline stages
Count = P_stage;        % The last stage contributing to Out

LPF_status = Pipeline_Type{3};
LPF_inv_status = Pipeline_Type{4};
Compensation = Pipeline_Type{2};
CLPTS = Pipeline_Type{1};
%LPF_status = 'OFF';       % Whether a LPF is used between the stages.
% LPF_inv_status = 'OFF';   % Whether the inverse of LPF is taken into account for digital output calculation.
%PZ_scaling = 'OFF';       % Scaling the poles and zeros of LPF between the pipeleine stages.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ***** S&H *******************

OSR = SH_feature(1);            % Over Sample Ratio

% *****************************
% ***** Amplifier *************

K = Amplifier_feature{1};              % Gain

% *****************************
% ***** ADC *******************

ADC_res = ADC_feature{1};      % ADC resolution
%ADC2_res = 10;
%ADC_res = [6;6];
%Max_INL = ADC_feature(2);        % Maximum INL (in LSB)
%Noise = 1;          % STD of noise added to quantizer levels

% *****************************
% ***** DAC *******************

DAC_res = DAC_feature{1};       % DAC resolution
Max_INL = DAC_feature{2};        % Maximum INL (in LSB)
%Noise = 1;          % STD of noise added to quantizer levels

% *****************************
% ***** Saturation *************

Sat_lim = SAT_ana;        % Saturation limit (dual polarity)
SAT_Dig = 2.^(ADC_res - 1);

%******************************
% ***** Predictor *************

% Everything is about the function 1/(1+kH(z)z^(-1))
%DC_gain = 1/50;        % DC gain
for i = 1:P_stage
    Order(i) = length(H_N{i}) - 1;          % order of filter
end
Order_LPF = length(LPF_N) - 1;  % Order of the LPF filters used between the pipeline stages
for i = 1:P_stage - 1
    try
Order_FF(i) = length(FF_N{i}) - 1;
    end
end
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
    X = Amp*sin(2*pi*fin*t); % Input Signal
elseif strcmp(In_type,'Chirp')
    X = Amp*chirp(t,0,P,fin); % Input Signal
elseif strcmp(In_type,'Multitone')
    X = sin(2*pi*fin*t) + sin(2*pi*fin/2*t) + sin(2*pi*fin/4*t) + sin(2*pi*3*fin/4*t); % Input Signal
    X = Amp*X/max(X);
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
DAC_Err = zeros(P_stage,L);            % DAC error
Out = zeros(1,L);               % Output digital
Out_f = zeros(1,L);             % Filtered Out
Out_bf = zeros(P_stage,L);  % A vector storing the digital number Out so that applies it to the inverse filter
Err = zeros(1,L);               % Whole ADC Error
G = zeros(P_stage,L);           % Adaptive gain


for p = 1:P_stage
    if CLPTS(p)
        for i = 1:L
            
            %%%%% Subtractor  %%%%%%%%%
            %i;
            E1(p,i) = S2(p,i) - DAC(p,i);
            %S2(i);
            %S2(i) - DAC(i);
            %%%%% Gain stage %%%%%%%%%%
            
            
            %G(i) = (1 - exp(-0.0001*i))*K;
            E2(p,i) = K(p)*E1(p,i);
            %K*E1(i)
            
            
            %%%%% Saturation %%%%%%%%%%
            
            SAT(p,i) = E2(p,i);
            
            
                if (E2(p,i) > Sat_lim(p))
                    SAT(p,i) = Sat_lim(p);
                    
                elseif (E2(p,i) < -Sat_lim(p))
                    SAT(p,i) = -Sat_lim(p);
                    
                end
           
            
            %%%%% Pipelining: LPF between the stages %%%%%%%%%%
            
            
            if (p < P_stage) % && (i < L)
                
                if LPF_status(p)
                    
                    %%%%%%%%%%%%%%%%%%%  LPF between the stages    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if i < Order_FF(p) + 1
                        S2(p + 1,i) = 0;
                        %S2(p,i + 1) = 0;
                    else
                        
                        S2(p + 1,i) = S2(p + 1,i - Order_FF(p):i - 1)*(flip(-FF_D{p}(2:end)))' + ...
                            SAT(p,i - Order_FF(p):i)*(flip(FF_N{p}(1:end)))';
                        
                        %Pred(i) = Pred(i) + H_N(1)*(2*ADC(i) - ADC(i -1));
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    
                    S2(p + 1,i) = SAT(p,i);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% ADC %%%%%%%%%%%%%%%%%%
            
            ADC(p,i) = floor(2^(ADC_res(p) - 1)*SAT(p,i)/Sat_lim(p));
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Saturation block %%%%%
            
            SAT2(p,i) = Pred(p,i);
            
                
                if (Pred(p,i) > 2^(ADC_res(p) - 1))
                    SAT2(p,i) = 2^(ADC_res(p) - 1);
                elseif (Pred(p,i) < -2^(ADC_res(p) - 1))
                    SAT2(p,i) = -2^(ADC_res(p) - 1);
                end

            
            %%%% DAC %%%%%%%%%%%%%%%%%%
            
            Y1 = SAT2(p,i)/2^(ADC_res(p) - DAC_res(p));
            Y2 = floor(Y1) + 0.5 - Max_INL(p)*sin(pi*Y1/2^(DAC_res(p) - 1));
            Y3 = floor(Y1) + 0.5;
            
            if i < L
                
                DAC(p,i + 1)     = Sat_lim(p)*Y2/2^(DAC_res(p) - 1);
                DAC_i(p,i + 1)   = Sat_lim(p)*Y3/2^(DAC_res(p) - 1);
                DAC_Err(p,i + 1) = DAC_i(p,i + 1) - Y1;
                
            end
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        for i = 1:L
            
            
            ADC(p,i) = floor(2^(ADC_res(p) - 1)*S2(p,i)/Sat_lim(p));
            %floor(2^(ADC_res(p) - 1)*S2(p,i)/Sat_lim(p))
            
            %%%% DAC %%%%%%%%%%%%%%%%%%
            
            Y1 = ADC(p,i)/2^(ADC_res(p) - DAC_res(p));
            Y2 = floor(Y1) + 0.5 - Max_INL(p)*sin(pi*Y1/2^(DAC_res(p) - 1));
            Y3 = floor(Y1) + 0.5;
                        
            DAC(p,i)   = Sat_lim(p)*Y2/2^(DAC_res(p) - 1);
            %Sat_lim(p)*Y2/2^(DAC_res(p) - 1)
            DAC_i(p,i) = Sat_lim(p)*Y3/2^(DAC_res(p) - 1);
                                 
            E1(p,i) = S2(p,i) - DAC(p,i);
            E2(p,i) = K(p)*E1(p,i);
            
            %%%%% Saturation %%%%%%%%%%
            
            SAT(p,i) = E2(p,i);
            
            
            if (E2(p,i) > Sat_lim(p))
                SAT(p,i) = Sat_lim(p);
                
            elseif (E2(p,i) < -Sat_lim(p))
                SAT(p,i) = -Sat_lim(p);
                
            end
            
            
            %%%%% Pipelining: LPF between the stages %%%%%%%%%%
            
            
            if p < P_stage % && i < L
                
                if LPF_status(p)
                    
                    %%%%%%%%%%%%%%%%%%%  LPF between the stages    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if i < Order_FF(p) + 1
                        S2(p + 1,i) = 0;
                        %S2(p,i + 1) = 0;
                    else
                        
                        S2(p + 1,i) = S2(p + 1,i - Order_FF(p):i - 1)*(flip(-FF_D{p}(2:end)))' + ...
                            SAT(p,i - Order_FF(p):i)*(flip(FF_N{p}(1:end)))';
                        
                        %Pred(i) = Pred(i) + H_N(1)*(2*ADC(i) - ADC(i -1));
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                else
                    %S2(p,i + 1) = SAT(p - 1,i);
                    S2(p + 1,i) = SAT(p,i);
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Pipeline Type 2 %%%%%%%%%%%%%%%%%%%%%%%%


% In this topology, the first stage is the proposed ADC (Closed-loop
% Two-step) and the second stage is a conventinal ADC.

%if Pipeline_Type == 2

% ADC2 = floor(2^(ADC2_res - 1)*SAT(1,:)/Sat_lim);
% %ADC2 = E2(1,:);

Out = Sat_lim(Count)*(ADC(Count,:) + 0.5)/(2^(ADC_res(Count) - 1));
if ~CLPTS(Count)
    Count = Count - 1;
    if P_stage > 1
    Out_bf(P_stage - 1,:) = Out;
    %Out_bf(P_stage - 1,:) = S2(P_stage,:);
    end
end

%Out = (ADC2)/K(1) + DAC(1,:);


%%%%%%%%%%%%%%%%%%  Digital domain processing %%%%%%%%%%%%%%%%%

%     if i - (Count -1) <1
%         Err(i) = 0;
%         Out(i) = 0;
%     else

%elseif Pipeline_Type == 1

%Out = (ADC(Count,:) + 0.5)/(2^(ADC_res(Count) - 1));

for pp = 1:Count
    p = Count - pp + 1;

        %%%%%%%%%%%%%%% Inverse Filter %%%%%%%%%%%%%%%%%%
    if  p < P_stage && LPF_inv_status(p)

            
            for i = 1:L
                
                if i < Order_FF(p) + 1
                    Out(i) = SAT(p,i);
                else
                    
                    Out(i) = Out(i - Order_FF(p):i - 1)*(flip(-IF_D{p}(2:end)))' + ...
                        Out_bf(p,i - Order_FF(p):i)*(flip(IF_N{p}(1:end)))';
                    
                end
                
                %                     S_O(i+1) = S2(2,i) - Out_bf(2,i);
                %                     if i < Order_LPF(p - 1) + 1
                %                         E2_est(i) = 0;
                %                     else
                %
                %                         E2_est(i) = E2_est(i - Order_LPF(p - 1):i - 1)*(flip(-LPF_inv_D{p - 1}(2:end)))' + ...
                %                             S_O(i - Order_LPF(p - 1)+1:i+1)*(flip(LPF_inv_N{p - 1}(1:end)))';
                %
                %                     end
            end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Compensation  %%%%%%%%%%%%%%%%%%
    if  p < P_stage && Compensation(p) && CLPTS(p)

            
            for i = 1:L
                
                if i < Order_FF(p) + 1
                    Out(i) = SAT(p,i);
                else
                    
                    Out(i) = Out(i - Order_FF(p):i - 1)*(flip(-IF_D{p}(2:end)))' + ...
                        Out_bf(p,i - Order_FF(p):i)*(flip(IF_N{p}(1:end)))';
                    
                end

            end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%% Inverse Filter %%%%%%%%%%%%%%%%%%
%     if  p > 1 && LPF_inv_status(p - 1)
% 
%             
%             for i = 1:L-1
%                 
%                 if i < Order_FF(p - 1) + 1
%                     Out(i) = SAT(p - 1,i);
%                 else
%                     
%                     Out(i) = Out(i - Order_FF(p - 1):i - 1)*(flip(-IF_D{p - 1}(2:end)))' + ...
%                         Out_bf(p,i - Order_FF(p - 1) + 1:i + 1)*(flip(IF_N{p - 1}(1:end)))';
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
%         
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t(10:10:end),Out,t(10:10:end),SAT(1,:))
legend('out','SAT1')

    Out = Out/K(p) + DAC_i(p,:);
    %Out = Out/K(p) + [ones(1,pp) DAC_i(p,1:end - pp)];
    %(p - Count)
    % + (p - Count)
    if p > 1
    Out_bf(p - 1,:) = Out;
    end
end
try
figure
plot(t(10:10:end),SAT(1,:),t(10:10:end),S2(2,:))
legend('SAT1','Input2')
end

%Out(i) = (ADC(2,i) + 0.5)/(2^(ADC_res(2) - 1))/K(1)/K(2) + DAC(2,i)/K(1) + DAC(1,i-1);

%         if i - (Count -1) <1
%             Err(i) = 0;
%             %Out(i) = 0;
%         else
%             Err(i) = Out(i) - S2(1,i-Count+1);
%         end
%Err = Out - [ones(1,P_stage - 1) S2(1,1:end - P_stage + 1)];
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
    if i < Order_LPF + 1
        Out_f(i + 1) = 0;
    else
        
        Out_f(i + 1) = Out_f(i - Order_LPF + 1:i)*(flip(-LPF_D(2:end)))' + ...
            Out(i - Order_LPF:i)*(flip(LPF_N(1:end)))';
        
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
In_status = zeros(1,P_stage);
for p = 1:P_stage
    Sat_status(p) = max(abs(E2(p,10*N_fft*OSR:L)));
    In_status(p)  = max(abs(S2(p,10*N_fft*OSR:L)));
end
%norm(Err(10*N_fft*OSR:L))

SNDR_t = 20*log10(norm(S2(1,10*N_fft*OSR:L))/norm(Err(10*N_fft*OSR:L)));
%     SNR2 = 20*log10(norm(S2(2,10*N_fft*OSR:L))/norm(Err2(10*N_fft*OSR:L)));
%norm(S2(N_fft*OSR:L))
%norm(Err(N_fft*OSR:L))

ENOB_t = (SNDR_t - 1.76)/6.02;
%     ENOB2(j) = (SNR2 - 1.76)/6.02;
% [ENOB,Amp]
%     [ENOB2(j),A(j)]
%ENOB2 = ENOB + log2(Sat_lim/A);





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
if strcmp(Plot_status{9,1},'Only Out')
    
    try
        figure
        plot(t(10:10:end),Out)
        title('Time domain output')
        xlabel('Time')
        %legend('Output Digital')
    end
elseif strcmp(Plot_status{9,1},'Out and Input')
    try
        figure
        plot(t(10:10:end),Out,t(10:10:end),X_S)
        title('Time domain input and output')
        xlabel('Time')
        legend('Output Digital','Input Signal')
    end
end
%
if strcmp(Plot_status{11,1},'Only Out')
    try
        figure
        plot(t(10:10:end),Out_f)
        title('Time domain filtered output')
        xlabel('Time')
        %legend('Filtered Output')
    end
elseif strcmp(Plot_status{11,1},'Out and Input')
    try
        figure
        plot(t(10:10:end),Out_f,t(10:10:end),X_S)
        title('Time domain filtered output and Input')
        xlabel('Time')
        legend('Filtered Output','Input Signal')
    end
end
%
if Plot_status{10,1}
    try
        figure
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(Out(L/2 + 1:L)))/max(abs(fft(Out(L/2 + 1:L))))))
        ylabel('dB')
        xlabel('Frequency')
        title('Digital Out FFT')
        %     L_Out_f = log2(length(Out_f))
    end
end
%
if Plot_status{12,1}
    try
        figure
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(Out_f(L/2 + 1:L)))/max(abs(fft(Out_f(L/2 + 1:L))))))
        ylabel('dB')
        xlabel('Frequency')
        title('Digital Out Filtered FFT')
        %     L_Out_f = log2(length(Out_f))
    end
end
%
for p = 1:P_stage
if Plot_status{13,p}
    try
        figure
        semilogy(abs(fft(E2(p,L/2 + 1:L))))
        title(['Amplifier output frequency domain: Stage',sprintf('%d',p)])
        xlabel('Frequency')
        legend('Amplifier FFT')
        %     L_E2 = log2(length(E2(1,:)))
    end
end
end
%
for p = 1:P_stage
if Plot_status{4,p}
    try
        figure
        plot(t(10:10:end),E2(p,:))
        title(['Amplifier output time domain',sprintf('%d',p)])
        xlabel('Time')
        legend('Amplifier output')
    end
end
end
%
for i = 1:P_stage
if Plot_status{1,p}
    try
        figure
        plot(t(10:10:end),Pred(p,:))
        xlabel('Time')
        title(['Prediction',sprintf('%d',p)])
    end
end
end
for i = 1:P_stage
if Plot_status{2,1}
    try
        figure
        plot(t(10:10:end),S2(p,:))
        xlabel('Time')
        title(['Input: Stage ',sprintf('%d',p)])
    end
end
end
%
if Plot_status{2,1} && Plot_status{3,1}
    try
        figure
        plot(t,X,t,S1)
        xlabel('Time')
        title('Input ans Sampled Signals')
        legend('Input', 'Sampled')
    end
elseif Plot_status{2,1}
    try
        figure
        plot(t,X)
        xlabel('Time')
        title('Input Signals')
    end
elseif Plot_status{3,1}
    try
        figure
        plot(t,S1)
        xlabel('Time')
        title('Sampled Signals')
    end
end
%
for p = 1:P_stage
if Plot_status{5,p}
    try
        figure
        plot(t(10:10:end),SAT(p,:))
        xlabel('Time')
        title(['Analog Saturation block',sprintf('%d',p)])
    end
end
end
%
for p = 1:P_stage
    try
if Plot_status{6,p}
    
        figure
        plot(t(10:10:end),SAT2(p,:))
        xlabel('Time')
        title(['Digital Saturation block',sprintf('%d',p)])
    
end
    end
end
%
for p = 1:P_stage
if Plot_status{7,p}
    try
        figure
        plot(t(10:10:end),ADC(p,:))
        xlabel('Time')
        ylabel('Digital Code')
        title(['ADC output',sprintf('%d',p)])
    end
end
end
%
for p = 1:P_stage
if Plot_status{8,p}
    try
        figure
        plot(t(10:10:end),DAC(p,:))
        xlabel('Time')
        title(['DAC output',sprintf('%d',p)])
    end
end
end
%

if Plot_status{14,1}
    for p = 1:P_stage
    try
        figure
        plot(t(10:10:end),S2(p,:),t(10:10:end),DAC(p,:))
        xlabel('Time')
        title(['Prediction vs actual input',sprintf('%d',p)])
        legend('Actual input','Prediction')
    end
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
try
    figure
    plot(Sat_status)
    xlabel('Stage number')
    title('Max of Amplifier output')
end
try
    figure
    plot(In_status)
    xlabel('Stage number')
    title('Max of Input')
end
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
FFT_D_O = abs(fft(Out(L/2 + 1:L)));
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

FFT_D_O = abs(fft(Out_f(L/2 + 1:L)));
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
end
