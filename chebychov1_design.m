%% ================================
%  3rd-Order Chebyshev Type I HPF
%  ================================
clear; clc;

%% Parameters
fs = 600;          % Sampling rate (Hz)
fc = 50;            % Desired digital cutoff (Hz)
n  = 3;             % Filter order
Rp = 0.5;           % Passband ripple (dB), typical: 0.2 to 1 dB
Nfft = 4096;        % For smooth freq response

Ts = 1/fs;

%% ===== Route A: Direct digital (cheby1 does BLT internally) =====
Wn = fc/(fs/2);                     % Normalized cutoff (0..1)
[bd, ad] = cheby1(n, Rp, Wn, 'high');  % Digital HP Chebyshev Type I

% Pole–zero map (digital)
figure;
zplane(bd, ad);
title(sprintf('Chebyshev-I HPF (Direct)  n=%d, Rp=%.2f dB, fc=%g Hz, fs=%g Hz', n, Rp, fc, fs));
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
disp('=== Digital 3rd-Order Chebyshev-I HPF (Direct) ===');
sys_dir

[zd, pd, kd] = tf2zpk(bd, ad);
disp('Zeros (Direct):'); disp(zd.');
disp('Poles (Direct):'); disp(pd.');
disp('Gain  (Direct):'); disp(kd);

num_str = poly2str(bd, 'z^{-1}');
den_str = poly2str(ad, 'z^{-1}');
fprintf('H_dir(z) = (%s) / (%s)\n\n', num_str, den_str);

%% ===== Route B: Explicit analog prototype + prewarp + BLT =====
% 1) Analog low-pass Chebyshev Type I prototype at 1 rad/s with Rp dB ripple
[za, pa, ka] = cheb1ap(n, Rp);      % analog LP prototype (wc=1 rad/s)
[ba_lp, aa_lp] = zp2tf(za, pa, ka);

% 2) Prewarp digital cutoff to analog (rad/s)
Omega_c = 2*fs*tan(pi*fc/fs);       % exact prewarp for BLT

% 3) Low-pass -> High-pass transform in s-domain
[ba_hp, aa_hp] = lp2hp(ba_lp, aa_lp, Omega_c);

% 4) Bilinear transform (s -> z)
[bb, ab] = bilinear(ba_hp, aa_hp, fs);

% Pole–zero map (digital, explicit BLT)
figure;
zplane(bb, ab);
title(sprintf('Chebyshev-I HPF (Analog prototype + BLT)  n=%d, Rp=%.2f dB', n, Rp));
grid on;

% Frequency response (linear) and comparison
[H2, f2] = freqz(bb, ab, Nfft, fs);

figure; hold on; grid on;
plot(f,  abs(H),  'LineWidth', 1.6, 'DisplayName', 'Direct (cheby1)');
plot(f2, abs(H2), 'LineWidth', 1.2, 'LineStyle','--', 'DisplayName', 'Analog+BLT');
xlabel('Frequency (Hz)'); ylabel('Magnitude (linear)');
title('Magnitude Response Comparison (linear scale)');
legend; xlim([0 fs/2]);
ylim([0 1.1*max([abs(H(:)); abs(H2(:))])]);

% Pretty-print H(z) in z^-1 form (Analog+BLT)
sys_blt = tf(bb, ab, Ts, 'variable', 'z^-1');
sys_blt = minreal(sys_blt, 1e-12);
disp('=== Digital 3rd-Order Chebyshev-I HPF (Analog + BLT) ===');
sys_blt

[zb, pb, kb] = tf2zpk(bb, ab);
disp('Zeros (BLT):'); disp(zb.');
disp('Poles (BLT):'); disp(pb.');
disp('Gain  (BLT):'); disp(kb);

num_str2 = poly2str(bb, 'z^{-1}');
den_str2 = poly2str(ab, 'z^{-1}');
fprintf('H_blt(z) = (%s) / (%s)\n\n', num_str2, den_str2);

%% (Optional) Spec-based design:
% If you have passband/stopband specs, solve Wn via cheb1ord first:
% Rp  = 0.5;     % passband ripple (dB)
% Rs  = 40;      % stopband attenuation (dB)
% fp  = 80;      % passband edge (Hz)
% fsb = 60;      % stopband edge (Hz)
% [n_min, Wn_min] = cheb1ord(fp/(fs/2), fsb/(fs/2), Rp, Rs, 'high');
% [b_opt, a_opt]  = cheby1(n_min, Rp, Wn_min, 'high');
