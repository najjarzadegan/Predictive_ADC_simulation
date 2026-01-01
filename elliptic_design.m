%% ============================================
%  3rd-Order Elliptic (Cauer) High-Pass Filter
%  ============================================
clear; clc;close all

%% Parameters
fs   = 1000;        % Sampling rate (Hz)
osr  = 3
fc   = fs/2/osr;          % Desired digital cutoff (Hz)
n    = 3;           % Filter order
Rp   = 0.5;         % Passband ripple (dB)
Rs   = 40;          % Stopband attenuation (dB)
Nfft = 4096;

Ts = 1/fs;

%% ===== Route A: Direct digital (ellip does BLT internally) =====
Wn = fc/(fs/2);                          % Normalized cutoff (0..1)
[bd, ad] = ellip(n, Rp, Rs, Wn, 'high'); % Digital HP Elliptic

% Pole–zero map (digital)
figure;
zplane(bd, ad);
title(sprintf('Elliptic HPF (Direct)  n=%d, Rp=%.2f dB, Rs=%g dB, fc=%g Hz, fs=%g Hz', ...
              n, Rp, Rs, fc, fs));
grid on;

% Frequency response (linear magnitude)
[H, f] = freqz(bd, ad, Nfft, fs);
figure;
plot(f, abs(H), 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (linear)');
title('Magnitude Response (Direct digital, linear scale)');
grid on; xlim([0 fs/2]);
ylim([0 1.1*max(abs(H))]);

% Phase (optional)
figure;
plot(f, unwrap(angle(H))*180/pi, 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Phase Response (Direct digital)');
grid on; xlim([0 fs/2]);

%% Pretty-print H(z) in z^-1 form (Direct)
sys_dir = tf(bd, ad, Ts, 'variable', 'z^-1');
sys_dir = minreal(sys_dir, 1e-12);
disp('=== Digital 3rd-Order Elliptic HPF (Direct) ===');
sys_dir

[zd, pd, kd] = tf2zpk(bd, ad);
disp('Zeros (Direct):'); disp(zd.');
disp('Poles (Direct):'); disp(pd.');
disp('Gain  (Direct):'); disp(kd);

num_str = poly2str(bd, 'z^{-1}');
den_str = poly2str(ad, 'z^{-1}');
fprintf('H_dir(z) = (%s) / (%s)\n\n', num_str, den_str);

%% ===== Route B: Explicit analog prototype + prewarp + BLT =====
% 1) Analog low-pass Elliptic prototype (wc = 1 rad/s, Rp dB pass ripple, Rs dB stop)
[za, pa, ka] = ellipap(n, Rp, Rs);   % zeros on jΩ, finite stopband zeros present
[ba_lp, aa_lp] = zp2tf(za, pa, ka);

% 2) Prewarp digital cutoff to analog (rad/s)
Omega_c = 2*fs*tan(pi*fc/fs);

% 3) Low-pass -> High-pass transform in s-domain
[ba_hp, aa_hp] = lp2hp(ba_lp, aa_lp, Omega_c);

% 4) Bilinear transform (s -> z)
[bb, ab] = bilinear(ba_hp, aa_hp, fs);

% Pole–zero map (digital, explicit BLT)
figure;
zplane(bb, ab);
title(sprintf('Elliptic HPF (Analog prototype + BLT)  n=%d, Rp=%.2f dB, Rs=%g dB', n, Rp, Rs));
grid on;

% Frequency response (linear) and comparison
[H2, f2] = freqz(bb, ab, Nfft, fs);

figure; hold on; grid on;
plot(f,  abs(H),  'LineWidth', 1.6, 'DisplayName', 'Direct (ellip)');
plot(f2, abs(H2), 'LineWidth', 1.2, 'LineStyle','--', 'DisplayName', 'Analog+BLT');
xlabel('Frequency (Hz)'); ylabel('Magnitude (linear)');
title('Magnitude Response Comparison (linear scale)');
legend; xlim([0 fs/2]);
ylim([0 1.1*max([abs(H(:)); abs(H2(:))])]);

% Pretty-print H(z) in z^-1 form (Analog+BLT)
sys_blt = tf(bb, ab, Ts, 'variable', 'z^-1');
sys_blt = minreal(sys_blt, 1e-12);
disp('=== Digital 3rd-Order Elliptic HPF (Analog + BLT) ===');
sys_blt

[zb, pb, kb] = tf2zpk(bb, ab);
disp('Zeros (BLT):'); disp(zb.');
disp('Poles (BLT):'); disp(pb.');
disp('Gain  (BLT):'); disp(kb);

num_str2 = poly2str(bb, 'z^{-1}');
den_str2 = poly2str(ab, 'z^{-1}');
fprintf('H_blt(z) = (%s) / (%s)\n\n', num_str2, den_str2);

%% (Optional) Spec-based minimized order:
% If you have pass/stop edges, solve order and Wn via ellipord:
% Rp  = 0.5;     % passband ripple (dB)
% Rs  = 40;      % stopband attenuation (dB)
% fp  = 80;      % passband edge (Hz)
% fsb = 60;      % stopband edge (Hz)
% [n_min, Wn_min] = ellipord(fp/(fs/2), fsb/(fs/2), Rp, Rs, 'high');
% [b_opt, a_opt]  = ellip(n_min, Rp, Rs, Wn_min, 'high');
