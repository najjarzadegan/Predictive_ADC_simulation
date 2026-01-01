%% ==============================================
%  Build custom H(z) from zeros/poles, print:
%    1) H(z)
%    2) G(z) = (1/H(z) - 1)*z
%    3) y[n] in terms of x[.] and y[.]
% ===============================================

clear; clc;


%% ---- Inputs (custom zeros/poles) ----
z0      = 0.65;             % real zero
z1      = 0.96;             % radius of complex-zero pair
theta_z = 30*pi/180;       % angle (radians) of complex-zero pair

p0      = -0.55;             % real pole
p1      = 0.9;             % radius of complex-pole pair
theta_p = 75*pi/180;       % angle (radians) of complex-pole pair

Ts = 1;  % sample time (discrete-time TF); change if you need a specific Ts
k=11.56
%% ---- Build zeros/poles and polynomials in z ----
z_zeros = [ z0, z1*exp(1j*theta_z), z1*exp(-1j*theta_z) ];
z_poles = [ p0, p1*exp(1j*theta_p), p1*exp(-1j*theta_p) ];

% b(z), a(z) with descending powers of z (variable = 'z')
b = poly(z_zeros);
a = poly(z_poles);

% H(z) in the 'z' variable
Hz = tf(b, a, Ts, 'variable', 'z');
Hz = minreal(Hz, 1e-12);

disp('=== H(z) from custom zeros/poles ===');
disp(Hz);
fprintf('H(z) = (%s) / (%s)\n\n', poly2str(b,'z'), poly2str(a,'z'));

%% ---- Compute G(z) = (1/H(z) - 1) * z ----
Gz = (1/k)*(1/Hz - 1) * tf([1 0], 1, Ts, 'variable', 'z');   % [1 0] is "z" in 'z' variable
Gz = minreal(Gz, 1e-12);

disp('=== G(z) = (1/H(z) - 1) * z ===');
disp(Gz);
[numG, denG] = tfdata(Gz, 'v');   % if you also want the raw polynomials:
fprintf('G(z) = (%s) / (%s)\n\n', poly2str(numG,'z'), poly2str(denG,'z'));

%% ---- Difference equation y[n] = ... (use z^{-1} form) ----
% Convert H(z) -> H(q) with q = z^{-1}:
% If H(z) = (b0 z^N + ... + bN)/(a0 z^N + ... + aN),
% then H(q) = (bN + b_{N-1} q + ... + b0 q^N)/(aN + ... + a0 q^N).
bz1 = fliplr(b);
az1 = fliplr(a);

% Normalize so a0 (in q-domain) = 1
bz1 = bz1 / az1(1);
az1 = az1 / az1(1);

% Now H(q) = (b0 + b1 q + ... + bN q^N) / (1 + a1 q + ... + aN q^N)
% Difference equation:
% y[n] = -a1 y[n-1] - a2 y[n-2] - ... - aN y[n-N] + b0 x[n] + b1 x[n-1] + ... + bN x[n-N]
N = max(numel(bz1), numel(az1)) - 1;
% pad (defensive)
bz1 = [bz1, zeros(1, N+1-numel(bz1))];
az1 = [az1, zeros(1, N+1-numel(az1))];

disp('=== Difference equation coefficients (z^{-1} form) ===');
fprintf('Numerator (b): '); fprintf('% .8g ', bz1); fprintf('\n');
fprintf('Denominator(a): '); fprintf('% .8g ', az1); fprintf('\n\n');

% Pretty print the equation
fprintf('y[n] = ');
for k = 1:N
    c = -az1(k+1);
    if abs(c) > 0
        fprintf('%+ .8g*y[n-%d] ', c, k);
    end
end
for k = 0:N
    c =  bz1(k+1);
    if abs(c) > 0
        if k == 0
            fprintf('%+ .8g*x[n] ', c);
        else
            fprintf('%+ .8g*x[n-%d] ', c, k);
        end
    end
end
fprintf('\n');

%% ---- (Optional) quick sanity plots ----
% Pole-zero map and magnitude response (linear)
figure; hold on; grid on; axis equal;
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k--','LineWidth',1,'DisplayName','Unit circle');

zH = roots(b); pH = roots(a);
plot(real(zH), imag(zH), 'bo', 'DisplayName','H zeros');
plot(real(pH), imag(pH), 'bx', 'DisplayName','H poles');
xlabel('Real'); ylabel('Imag'); title('PZ map of H(z)'); legend; hold off;

figure; 
[Hf,w] = freqz(fliplr(b), fliplr(a), 2048);   % freqz uses z^{-1} form; flip b,a
plot(w, abs(Hf), 'LineWidth',1.5);
xlabel('\omega (rad/sample)'); ylabel('|H(e^{j\omega})| (linear)');
title('Magnitude response of H(z)'); grid on;
