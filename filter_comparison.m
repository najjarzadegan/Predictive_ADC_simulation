%% ============================================================
%  High-Pass (HP) Filter Comparison, Order 3, Common OSR cutoff
%  Butterworth, Chebyshev I, Chebyshev II, Elliptic, and Custom
%  - Cutoff in rad/sample: wc = pi/OSR  (fs/2 <-> pi)
%  - z^-1 form printed with constant terms normalized to 1
%  - Common PZ map and magnitude response (linear)
%  - Metric computed in rad/sample with HP bands:
%       Passband: [wc, pi]
%       Stopband: [0, wc]
% ============================================================

clear; clc;close all

%% ---------------- User Parameters ----------------
fs   = 1000;      % Sampling rate [Hz] (used only for pretty TF printing)
OSR1  = 2;         % Oversampling ratio  -> wc = pi/OSR
OSR2  = 5;
n    = 3;         % Filter order (3)

% Chebyshev / Elliptic specs
Rp   = 0.5;       % passband ripple [dB]   (Cheby1/Ellip)
Rs   = 40;        % stopband atten [dB]    (Cheby2/Ellip)

% Custom filter (digital z-plane)
z0      = 0.65;              % real zero
z1      = 0.96;              % radius of complex-zero pair
theta_z = 30*pi/180;        % angle (radians) of complex-zero pair

p0      = -0.55;              % real pole
p1      = 0.9;              % radius of complex-pole pair
theta_p = 75*pi/180;        % angle (radians) of complex-pole pair

% p0      = 0.6;              % real zero
% p1      = 0.9;              % radius of complex-zero pair
% theta_p = 30*pi/180;        % angle (radians) of complex-zero pair

% Metric variables
M = 7;      % user-chosen
N = 6;      % user-chosen
K = 11.56;      % user-chosen
A_p = 1.0;  % user-chosen

% Frequency grid for evaluation/plots (rad/sample)
Nfft = 8192;    % dense for smooth curves

%% ---------------- Derived Parameters ----------------
Ts = 1/fs;
wc = pi/OSR1;              % HP cutoff in rad/sample
wc2 = pi/OSR2;
Wn = wc/pi;               % normalized cutoff for design (0..1) for MATLAB design fns

%% ---------------- Helper functions ----------------
% Normalize so constant terms in z^-1 form are 1 (independently)
normalize_const_ones = @(b,a) deal( b./b(1), a./a(1) );

% Pretty TF print


% Mag response on [0, pi] (rad/sample)
mag_resp = @(b,a) freqz(b, a, Nfft);   % returns H, w (rad/sample)

%% ---------------- Design standard HIGH-PASS filters ----------------
% Butterworth
[b_but, a_but] = butter(n, Wn, 'high');

% Chebyshev Type I
[b_c1,  a_c1 ] = cheby1(n, Rp, Wn, 'high');

% Chebyshev Type II
[b_c2,  a_c2 ] = cheby2(n, Rs, Wn, 'high');

% Elliptic (Cauer)
[b_ell, a_ell] = ellip(n, Rp, Rs, Wn, 'high');

%% ---------------- Custom filter from specified z/p ----------------
zeros_custom = [ z0, z1*exp(1j*theta_z), z1*exp(-1j*theta_z) ];
poles_custom = [ p0, p1*exp(1j*theta_p), p1*exp(-1j*theta_p) ];
b_cus = poly(zeros_custom);
a_cus = poly(poles_custom);

%% ---------------- Normalize constant terms to 1 (z^-1 form) ----------------
[b_but, a_but] = normalize_const_ones(b_but, a_but);
[b_c1,  a_c1 ] = normalize_const_ones(b_c1,  a_c1 );
[b_c2,  a_c2 ] = normalize_const_ones(b_c2,  a_c2 );
[b_ell, a_ell] = normalize_const_ones(b_ell, a_ell);
[b_cus, a_cus] = normalize_const_ones(b_cus, a_cus);

%% ---------------- Print transfer functions ----------------
disp('Cutoff (rad/sample) from OSR:'); disp(wc);
disp('Note: Each TF is shown in z^{-1} form with constant terms normalized to 1 (num & den).');
disp(' ');

Ts = 1/fs;

print_tf('Butterworth HP (n=3)', b_but, a_but, Ts);
print_tf('Chebyshev-I HP (n=3, Rp)', b_c1, a_c1, Ts);
%print_tf('Chebyshev-II HP (n=3, Rs)', b_c2, a_c2, Ts);
print_tf('Elliptic HP (n=3, Rp, Rs)', b_ell, a_ell, Ts);
print_tf('Custom HP (user z/p)', b_cus, a_cus, Ts);

%% ---------------- Common Pole–Zero Map ----------------
figure; hold on; grid on; axis equal;
% Unit circle
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1, 'DisplayName','Unit circle');


plot_pz(b_but,a_but,'b','Butter');
plot_pz(b_c1 ,a_c1 ,'g','Cheb1');
plot_pz(b_c2 ,a_c2 ,'m','Cheb2');
plot_pz(b_ell,a_ell,'c','Ellip');
plot_pz(b_cus,a_cus,'r','Custom');

xlabel('Real'); ylabel('Imag');
title(sprintf('Pole–Zero Map (All HP, n=%d, \\omega_c=%.4f rad/sample)', n, wc));
legend('Location','eastoutside');
hold off;

%% ---------------- Magnitude Responses (linear, rad/sample) ----------------
[H_but, w]  = mag_resp(b_but,a_but);    % w in [0, pi] rad/sample
[H_c1,  ~]  = mag_resp(b_c1, a_c1 );
[H_c2,  ~]  = mag_resp(b_c2, a_c2 );
[H_ell, ~]  = mag_resp(b_ell,a_ell);
[H_cus, ~]  = mag_resp(b_cus,a_cus);

%% ---------------- Metric computation (rad/sample) ----------------
% HP bands:
%   Passband: [wc, pi]
%   Stopband: [0, wc]
pass_idx = (w >= wc2 - eps);
stop_idx = (w <= wc2 + eps);
all_indx = (w <= pi + eps);

% Helper: max magnitude in *passband*
pb_max = @(H) max(abs(H(stop_idx)));

% Helper: L2 norm over stopband of 1 - H, integrated over rad/sample
l2_stop = @(H) sqrt( trapz(w(all_indx), abs(1 - H(all_indx)).^2) / pi);

% Evaluate per filter
filters = {'Butter','Cheb1','Cheb2','Ellip','Custom'};
Hs      = { H_but,  H_c1,   H_c2,   H_ell,  H_cus  };

pbm = zeros(numel(filters),1);
l2n = zeros(numel(filters),1);
J   = zeros(numel(filters),1);

for i = 1:numel(filters)
    Hi = Hs{i};
    pbm(i) = pb_max(Hi);
    l2n(i) = l2_stop(Hi)
    disp(filters{i})
    term1 = 2^(-M)
    term2 = A_p * pbm(i)
    term3 = (1/K) * sqrt( 3 * ( 4^(-N) + K^2 * 4^(-M) ) ) * l2n(i)

    J(i) = term1 + term2 + term3;
end

%% ---------------- Plot responses with metric annotations ----------------
figure; hold on; grid on;
plot(w, abs(H_but), 'b', 'LineWidth', 1.6, 'DisplayName', sprintf('Butter  (J=%.6g)', J(1)));
plot(w, abs(H_c1 ), 'g', 'LineWidth', 1.6, 'DisplayName', sprintf('Cheb1   (J=%.6g)', J(2)));
plot(w, abs(H_c2 ), 'm', 'LineWidth', 1.6, 'DisplayName', sprintf('Cheb2   (J=%.6g)', J(3)));
plot(w, abs(H_ell), 'c', 'LineWidth', 1.6, 'DisplayName', sprintf('Ellip   (J=%.6g)', J(4)));
plot(w, abs(H_cus), 'r', 'LineWidth', 1.8, 'DisplayName', sprintf('Custom  (J=%.6g)', J(5)));

xline(wc2, 'k--', 'LineWidth', 1.2, 'DisplayName','\omega_c');

xlabel('\omega (rad/sample)'); ylabel('Magnitude (linear)');
title(sprintf('HP Magnitude Responses (linear)  OSR=%g, \\omega_c=\\pi/OSR=%.4f', OSR2, wc2));
xlim([0 pi]);
ylim([0 1.1*max([abs(H_but);abs(H_c1);abs(H_c2);abs(H_ell);abs(H_cus)])]);
legend('Location','southwest');

% Optional: add small text labels near end-of-band
x_annot = 0.96*pi;     % near Nyquist
yvals = [abs(H_but(end)), abs(H_c1(end)), abs(H_c2(end)), abs(H_ell(end)), abs(H_cus(end))];
dy = 0.05; y0 = min(yvals);
texts = {sprintf('J_Butter=%.6g',J(1)), sprintf('J_Cheb1=%.6g',J(2)), ...
         sprintf('J_Cheb2=%.6g',J(3)), sprintf('J_Ellip=%.6g',J(4)), ...
         sprintf('J_Custom=%.6g',J(5))};
for i=1:numel(texts)
    text(x_annot, min(1,max(0,y0 + (i-1)*dy)), texts{i}, 'FontSize', 10, 'Color', 'k');
end
hold off;

%% ---------------- Show metric results as a table ----------------
T = table(filters.', pbm, l2n, J, ...
    'VariableNames', {'Filter','PassbandMaxMag','L2_Stop_1minusH','MetricValue'});
disp('=== Metric results (HP bands; rad/sample) ===');
disp(T);

function print_tf(name, b, a, Ts)
    sys = tf(b, a, Ts, 'variable', 'z^-1');
    sys = minreal(sys, 1e-12);
    disp(['=== ' name ' ===']);
    disp(sys);
    fprintf('H(z) = (%s) / (%s)\n\n', poly2str(b,'z^{-1}'), poly2str(a,'z^{-1}'));
end

function plot_pz(b, a, clr, name)
    % Plot zeros
    z = roots(b);
    plot(real(z), imag(z), 'o', ...
        'Color', clr, 'MarkerSize', 6, 'LineWidth', 1.4, ...
        'DisplayName', [name ' zeros']);
    % Plot poles
    p = roots(a);
    plot(real(p), imag(p), 'x', ...
        'Color', clr, 'MarkerSize', 6, 'LineWidth', 1.4, ...
        'DisplayName', [name ' poles']);
end
