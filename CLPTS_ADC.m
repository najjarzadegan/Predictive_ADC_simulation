% Closed-Loop Predictive two step (CLPTS) ADC with shared ADC

function [Fs,BW,SNDR,ENOB,SNDR_f,ENOB_f,SNDR_t,ENOB_t,SAT_Dig] = CLPTS_ADC(Plot_status,...
          DAC_feature, ADC_feature, SAT_ana, Amplifier_feature,...
          SH_feature, Input_feature,In_Type,H_N,H_D,LPF_N,LPF_D,Out_type)

% Plot_status      : Determines which figures/signals are plotted.
% DAC_feature      : Denotes the DAC specifications.
% ADC_feature      : Denotes the ADC specifications.
% SAT_ana          : Denotes the Analog Saturation Limit.
% Amplifier_feature: Denotes the Amplifier specifications.
% SH_feature       : Denotes the Sample/Hold specifications.
% Input_feature    : Denotes the input signal specifications.
% In_Type          : Denotes the input type (sine, chirp, ...).
% H_N              : Predictor filters (Numirator).
% H_D              : Predictor filters (Denominator).
% LPF_N            : output LPF (Numirator).
% LPF_D            : output LPF (Denominator).
% Out_type         : How to generate the output

% Fs               : Sampling frequency
% BW               : Signal BandWidth
% SNDR             : Output SNDR (based on FFT)
% ENOB             : Output ENOB (based on FFT)
% SNDR_f           : Output SNDR after LPF (based on FFT)
% ENOB_f           : Output ENOB after LPF (based on FFT)
% SNDR_t           : Output SNDR before LPF (Time domain calculation)
% ENOB_t           : Output ENOB before LPF (Time domain calculation)
% SAT_Dig          : Denotes the Digital Saturation Limit.

tic                                 % Start time
clc                                 % Clean the workspace
close all                           % Close all the figures
%clear all                          % Clear all the variables (not
                                    % necessary)

% The system consists of a Sample and Hold (S&H), a subtractor, an ADC, a
% DAC and a digital filter working as predictor.

%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
% ***** Input signal **********

Amp = Input_feature{2};             % Amplitude
fin = Input_feature{3};             % frequency
In_type = In_Type;

% *****************************
% ***** Pipeline **************

P_stage = 1;                        % The number of pipeline stages
Count = P_stage;                    % The last stage contributing to Out

% *****************************
% ***** S&H *******************

OSR = SH_feature(1);                % Over Sample Ratio
N_fft = SH_feature(2);              % # of FFT points
R = SH_feature(3);                  % # of repearitions (R*M cycles) where M is the number of cycles we need for N_fft samples.

% *****************************
% ***** Amplifier *************

K = Amplifier_feature(1)*ones(P_stage,1);       % Gain

% *****************************
% ***** ADC *******************

ADC_res = ADC_feature(1)*ones(P_stage,1);       % ADC resolution
%Max_INL = ADC_feature(2);                      % Maximum INL (in LSB)
%Noise = 1;                                     % STD of the noise added to the quantizer levels

% *****************************
% ***** DAC *******************

DAC_res = DAC_feature(1)*ones(P_stage,1);       % DAC resolution
Max_INL = DAC_feature(2)*ones(P_stage,1);       % Maximum INL (in LSB)
%Noise = 1;                                     % STD of the noise added to the quantizer levels

% *****************************
% ***** Saturation ************

Sat_lim = SAT_ana*ones(P_stage,1);              % ADC Saturation limit (dual polarity)
SAT_Dig = 2.^(ADC_res - 1);                     % Predictor Saturation limit (dual polarity)

%******************************
% ***** Predictor *************

% Everything is about the function 1/(1+kH(z)z^(-1))

Order = length(H_N{1}) - 1;        % order of the filter
Order_LPF = length(LPF_N{1}) - 1;  % Order of the output LPF filter

% *****************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%% Filter Generation %%%%%%%%%%%%%%%%%%%%%%

% The predictor filter and the LPF are generated in another m-file and the
% transfer functions (H_N, H_D, LPF_N, LPF_D) are loaded to this function.

% *****************************

%%%%%%%%%%%% Parts %%%%%%%%%%%%%%%%%%%%

% ***** Input Signal *******************

M = N_fft/2 - 1;                % The number of cycles we need for N_fft samples.
Tin = 1/fin;                    % Input period.
P = Tin*M*R;                    % Simulation period.

Ts_Nyq = Tin*M/N_fft;           % Sampling Period for nyquist rate sampling
Ts_OS = Ts_Nyq/OSR ;            % Sampling Period for oversampling
Fs = 1/Ts_OS;                   % Samplig freq
BW = Fs/2/OSR;                  % Nyquist bandwidth

step = Ts_OS/10;                % Time step
t = step:step:P;                % Time vector (number of points is N_fft*R*OSR*10)

if strcmp(In_type,'Sine')
    X = Amp*sin(2*pi*fin*t);    % Input Signal = Sine
elseif strcmp(In_type,'Chirp')
    X = Amp*chirp(t,0,P,fin);   % Input Signal = Chirp
elseif strcmp(In_type,'Multitone')
    X = sin(2*pi*fin*t) + sin(2*pi*fin/2*t) + sin(2*pi*fin/4*t) + sin(2*pi*3*fin/4*t); 
    X = Amp*X/max(X);           % Input Signal = Multitone
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

E1 = zeros(P_stage,L);          % Loop Error
E2 = zeros(P_stage,L);          % Error after amplification
SAT = zeros(P_stage,L);         % Saturation block output
ADC = zeros(P_stage,L);         % ADC output
Pred = zeros(P_stage,L);        % Predictor output
SAT2 = zeros(P_stage,L);        % saturation output
DAC = zeros(P_stage,L);         % Initializing DAC
DAC_i = zeros(P_stage,L);       % Initializing DAC
Out = zeros(1,L);               % Output digital
Out_f = zeros(1,L);             % Filtered Out
Err = zeros(1,L);               % Whole ADC Error


for i = 1:L
    
    %%%%% Subtractor  %%%%%%%%%
    
    E1(:,i) = S2(:,i) - DAC(:,i);

    %%%%% Gain stage %%%%%%%%%%
    
    E2(:,i) = K.*E1(:,i);
    
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
   
    % Nonidealities like INL and noise should be added later
    
    
    %%%% Predictor %%%%%%%%%%%%
    
    for p = 1:P_stage
        

        if i < Order(p) + 1
            
            Pred(p,i) = 0;
            
        else
            
            Pred(p,i) = Pred(p,i - Order(p):i - 1)*(flip(-H_D{p}(2:end)))' + ...
                        ADC(p,i - Order(p) + 1:i)*(flip(H_N{p}(2:end)))';
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Saturation block (Digital) %%%%%
    
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
        
        DAC(:,i+1)   = Sat_lim.*Y2./2.^(DAC_res - 1);       % Actual DAC
        DAC_i(:,i+1) = Sat_lim.*Y3./2.^(DAC_res - 1);       % Ideal DAC
        
    end
      
end

%%%%%%%%%%%%  Digital domain processing %%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Out_type, 'Two-Step')
    
    Out = (ADC(1,:) + 0.5)/(2^(ADC_res(1) - 1))/K(1) + DAC_i(1,:);    % Constructing the output in a way a Two-step ADC does.
    
elseif strcmp(Out_type, 'Delta-Modulator')
    
    Out = ADC;                                                        % Constructing the output in a way a Delta-Modulator does.

else
    
    disp('Error in choosing the output construction method')
    
end

Err = Out - S2(1,:);                                                  % Conversion Error

%%%%%%%%%%%%%%%%%%%  Output Filter   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:L - 1
    
    if i < Order_LPF(P_stage) + 1
        
        Out_f(i + 1) = 0;       % Initial values for the output
        
    else
        
        Out_f(i + 1) = Out_f(i - Order_LPF(P_stage) + 1:i)*(flip(-LPF_D{P_stage}(2:end)))' + ...
            Out(i - Order_LPF(p):i)*(flip(LPF_N{P_stage}(1:end)))';
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sat_status = zeros(1,P_stage);


SNDR_t = 20*log10(norm(S2(1,10*N_fft*OSR:L))/norm(Err(10*N_fft*OSR:L)));

ENOB_t = (SNDR_t - 1.76)/6.02;

%
if strcmp(Plot_status{9},'Only Out')
    
    try
        figure
        plot(t(10:10:end),Out)
        title('Time domain output')
        xlabel('Time')
        %legend('Output Digital')
    end
    
elseif strcmp(Plot_status{9},'Out and Input')
    
    try
        figure
        plot(t(10:10:end),Out,t(10:10:end),X_S)
        title('Time domain input and output')
        xlabel('Time')
        legend('Output Digital','Input Signal')
    end
end
%
if strcmp(Plot_status{11},'Only Out')
    
    try
        figure
        plot(t(10:10:end),Out_f)
        title('Time domain filtered output')
        xlabel('Time')
        %legend('Filtered Output')
    end
    
elseif strcmp(Plot_status{11},'Out and Input')
    
    try
        figure
        plot(t(10:10:end),Out_f,t(10:10:end),X_S)
        title('Time domain filtered output and Input')
        xlabel('Time')
        legend('Filtered Output','Input Signal')
    end
end
%
if Plot_status{10}                          % Digital Out FFT
    
    try
        figure
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(Out(L/2 + 1:L)))/max(abs(fft(Out(L/2 + 1:L))))))
        ylabel('dB')
        xlabel('Frequency')
        title('Digital Out FFT')
        
    end
end
%
if Plot_status{12}                          % Digital Out Filtered FFT
    
    try
        figure
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(Out_f(L/2 + 1:L)))/max(abs(fft(Out_f(L/2 + 1:L))))))
        ylabel('dB')
        xlabel('Frequency')
        title('Digital Out Filtered FFT')
        
    end
end
%
if Plot_status{13}                          % Amplifier Output (Freq domain)
    
    try
        figure
        %semilogy(abs(fft(E2(1,L/2 + 1:L))))
        plot(2*Fs/L*(1:L/2),20*log10(abs(fft(E2(1,L/2 + 1:L)))/max(abs(fft(E2(1,L/2 + 1:L))))))
        title('Amplifier output frequency domain')
        ylabel('dB')
        xlabel('Frequency')
        legend('Amplifier FFT')
        
    end
end
%
if Plot_status{4}                           % Amplifier Output (Time domain)
    
    try
        figure
        plot(t(10:10:end),E2(1,:))
        title('Amplifier output time domain')
        xlabel('Time')
        legend('Amplifier output')
    end
end
%
if Plot_status{1}                           % Prediction
    
    try
        figure
        plot(t(10:10:end),Pred(1,:))
        xlabel('Time')
        title('Prediction')
    end
end
%
if Plot_status{2} && Plot_status{3}         % Comparing input and sampled
    
    try
        figure
        plot(t,X,t,S1)
        xlabel('Time')
        title('Input and Sampled Signals')
        legend('Input', 'Sampled')
    end
    
elseif Plot_status{2}                       % Input
    
    try
        figure
        plot(t,X)
        xlabel('Time')
        title('Input Signals')
    end
    
elseif Plot_status{3}                       % Sampled input
    
    try
        figure
        plot(t,S1)
        xlabel('Time')
        title('Sampled Signals')
    end
end
%
if Plot_status{5}                           % Analog saturation
    
    try
        figure
        plot(t(10:10:end),SAT(1,:))
        xlabel('Time')
        title('Analog Saturation block')
    end
end
%
if Plot_status{6}                           % Digital satiration
    
    try
        figure
        plot(t(10:10:end),SAT2(1,:))
        xlabel('Time')
        title('Digital Saturation block')
    end
end
%
if Plot_status{7}                           % ADC output
    
    try
        figure
        plot(t(10:10:end),ADC(1,:))
        xlabel('Time')
        ylabel('Digital Code')
        title('ADC output')
    end
end
%
if Plot_status{8}                           % DAC output
    
    try
        figure
        plot(t(10:10:end),DAC(1,:))
        xlabel('Time')
        title('DAC output')
    end
end
%
if Plot_status{14}                          % Comparing Input and prediction
    
    try
        figure
        plot(t(10:10:end),X_S,t(10:10:end),DAC)
        xlabel('Time')
        title('Prediction vs actual input')
        legend('Actual input','Prediction')
    end
end


%%%%%%%%% Calculating SNDR and ENOB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These calculations are close to the result from time-domain calculation
% only for single and milti tone sinusoidal signal. 

%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%

FFT_D_O = abs(fft(Out(L/2 + 1:L)));                 % FFT on the second half of the data

Max_bin = max(FFT_D_O(1:L/4));                      % Max bin = Input signal
Fm = find(FFT_D_O(1:L/4) > Max_bin/2);               % Input frequency

FFT_D_O_new = FFT_D_O(1:L/4);                       % 0 to Fs/2
FFT_D_O_new(Fm) = 0;                                % Removing the input

SNDR = 10*log10((norm(FFT_D_O(Fm)))^2/sum(FFT_D_O_new.^2));

ENOB = (SNDR - 1.76)/6.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FILTERED OUTPUT %%%%%%%%%%%%%%

FFT_D_O = abs(fft(Out_f(L/2 + 1:L)));

Max_bin = max(FFT_D_O(1:L/4));
Fm = find(FFT_D_O(1:L/4) > Max_bin/2);

FFT_D_O_new = FFT_D_O(1:L/4);
FFT_D_O_new(Fm) = 0;

SNDR_f = 10*log10((norm(FFT_D_O(Fm)))^2/sum(FFT_D_O_new.^2));

ENOB_f = (SNDR_f - 1.76)/6.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc                                                 % End of the simulation
