%% Parameters
fs = 1000;         % sampling rate (Hz)
fc = 60;           % desired digital cutoff (Hz)
n  = 3;            % filter order (3rd order)

%% ===== Route A: Direct digital (butter does BLT internally) =====
Wn = fc/(fs/2);                     % normalized cutoff (0..1)
[bd, ad] = butter(n, Wn, 'high');   % digital HP Butterworth

%% ---- Pretty-print H(z) in z^-1 form (Direct) ----
Ts = 1/fs;  % sample time

sys_dir = tf(bd, ad, Ts, 'variable', 'z^-1');   % discrete-time TF in z^-1
sys_dir = minreal(sys_dir, 1e-12);              % clean tiny coefficients
disp('=== Digital 3rd-Order Butterworth HP (Direct) ===');
sys_dir                                           % prints H(z)

% Also show zeros/poles/gain (z-domain)
[zd, pd, kd] = tf2zpk(bd, ad);
disp('Zeros (Direct):'); disp(zd.');
disp('Poles (Direct):'); disp(pd.');
disp('Gain  (Direct):'); disp(kd);

% Optional: print as polynomial strings (z^{-1} powers)
num_str = poly2str(bd, 'z^{-1}');
den_str = poly2str(ad, 'z^{-1}');
fprintf('H_dir(z) = (%s) / (%s)\n\n', num_str, den_str);

%% ---- Pretty-print H(z) in z^-1 form (Analog+BLT) ----
sys_blt = tf(bb, ab, Ts, 'variable', 'z^-1');
sys_blt = minreal(sys_blt, 1e-12);
disp('=== Digital 3rd-Order Butterworth HP (Analog + BLT) ===');
sys_blt

[zb, pb, kb] = tf2zpk(bb, ab);
disp('Zeros (BLT):'); disp(zb.');
disp('Poles (BLT):'); disp(pb.');
disp('Gain  (BLT):'); disp(kb);

num_str2 = poly2str(bb, 'z^{-1}');
den_str2 = poly2str(ab, 'z^{-1}');
fprintf('H_blt(z) = (%s) / (%s)\n', num_str2, den_str2);


% Pole–zero map (digital)
figure; 
zplane(bd, ad); 
title(sprintf('Digital 3rd-Order Butterworth HP (Direct)  fc=%g Hz, fs=%g Hz', fc, fs));
grid on;

% Frequency response (digital)
Nfft = 2048;
[H, f] = freqz(bd, ad, Nfft, fs);   % frequency vector in Hz
figure;
plot(f, abs(H), 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude Response (Direct digital)');
grid on; xlim([0 fs/2]);

figure;
plot(f, unwrap(angle(H))*180/pi, 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Phase Response (Direct digital)');
grid on; xlim([0 fs/2]);

%% ===== Route B: Explicit s-domain + prewarping + bilinear =====
% 1) Analog low-pass prototype at 1 rad/s
[za, pa, ka] = buttap(n);               % prototype Butterworth (LP, wc=1 rad/s)
[ba_lp, aa_lp] = zp2tf(za, pa, ka);

% 2) Prewarp the desired digital cutoff to analog rad/s
Omega_c = 2*fs*tan(pi*fc/fs);           % prewarped analog cutoff (rad/s)

% 3) Low-pass -> High-pass frequency transform in s-domain
[ba_hp, aa_hp] = lp2hp(ba_lp, aa_lp, Omega_c);

% 4) Bilinear transform (s -> z)
[bb, ab] = bilinear(ba_hp, aa_hp, fs);

%% ---- Pretty-print H(z) in z^-1 form (Direct) ----
Ts = 1/fs;  % sample time

sys_dir = tf(bd, ad, Ts, 'variable', 'z^-1');   % discrete-time TF in z^-1
sys_dir = minreal(sys_dir, 1e-12);              % clean tiny coefficients
disp('=== Digital 3rd-Order Butterworth HP (Direct) ===');
sys_dir                                           % prints H(z)

% Also show zeros/poles/gain (z-domain)
[zd, pd, kd] = tf2zpk(bd, ad);
disp('Zeros (Direct):'); disp(zd.');
disp('Poles (Direct):'); disp(pd.');
disp('Gain  (Direct):'); disp(kd);

% Optional: print as polynomial strings (z^{-1} powers)
num_str = poly2str(bd, 'z^{-1}');
den_str = poly2str(ad, 'z^{-1}');
fprintf('H_dir(z) = (%s) / (%s)\n\n', num_str, den_str);

%% ---- Pretty-print H(z) in z^-1 form (Analog+BLT) ----
sys_blt = tf(bb, ab, Ts, 'variable', 'z^-1');
sys_blt = minreal(sys_blt, 1e-12);
disp('=== Digital 3rd-Order Butterworth HP (Analog + BLT) ===');
sys_blt

[zb, pb, kb] = tf2zpk(bb, ab);
disp('Zeros (BLT):'); disp(zb.');
disp('Poles (BLT):'); disp(pb.');
disp('Gain  (BLT):'); disp(kb);

num_str2 = poly2str(bb, 'z^{-1}');
den_str2 = poly2str(ab, 'z^{-1}');
fprintf('H_blt(z) = (%s) / (%s)\n', num_str2, den_str2);

% Pole–zero map (digital, explicit BLT)
figure; 
zplane(bb, ab); 
title(sprintf('Digital 3rd-Order Butterworth HP (Analog + BLT)  fc=%g Hz, fs=%g Hz', fc, fs));
grid on;

% Frequency response (compare to Route A)
[H2, f2] = freqz(bb, ab, Nfft, fs);

figure; hold on; grid on;
plot(f,  abs(H),  'LineWidth', 1.5, 'DisplayName','Direct (butter)');
plot(f2, abs(H2), 'LineWidth', 1.2, 'LineStyle','--', 'DisplayName','Analog+BLT');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude Response Comparison');
legend; xlim([0 fs/2]);

%% Optional: show poles/zeros numerically
[z_dir, p_dir, k_dir] = tf2zpk(bd, ad);
disp('Direct digital (butter):');    disp(table(z_dir, p_dir));
[z_blt, p_blt, k_blt] = tf2zpk(bb, ab);
disp('Analog+BLT:');                 disp(table(z_blt, p_blt));
