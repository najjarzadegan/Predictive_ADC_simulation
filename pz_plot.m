% =====================================================
% Compare two transfer functions (H1 and H2)
% =====================================================

% ---- Parameters: First system ----
z0_1 = 0.65;
p0_1 = -0.55;
z1_1 = 0.96;
p1_1 = 0.9;
theta_z1 = 30*pi/180;
theta_p1 = 75*pi/180;

% ---- Parameters: Second system ----
z0_2 = 0.99;
p0_2 = 0.8;
z1_2 = 0.99;
p1_2 = 0.8;
theta_z2 = 36*pi/180;
theta_p2 = 36*pi/180;

% ---- Compute zeros and poles ----
zeros1 = [z0_1, z1_1*exp(-1j*theta_z1), z1_1*exp(1j*theta_z1)];
poles1 = [p0_1, p1_1*exp(-1j*theta_p1), p1_1*exp(1j*theta_p1)];

zeros2 = [z0_2, z1_2*exp(-1j*theta_z2), z1_2*exp(1j*theta_z2)];
poles2 = [p0_2, p1_2*exp(-1j*theta_p2), p1_2*exp(1j*theta_p2)];

% ---- Compute coefficients ----
b1 = poly(zeros1);
a1 = poly(poles1);
b2 = poly(zeros2);
a2 = poly(poles2);

% =====================================================
% Plot 1: Pole–Zero Map (custom color control)
% =====================================================
figure; hold on; grid on; axis equal;

% Plot unit circle
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);

% System 1: blue
plot(real(zeros1), imag(zeros1), 'bo', 'MarkerSize', 8, 'LineWidth', 1.8);
plot(real(poles1), imag(poles1), 'bx', 'MarkerSize', 8, 'LineWidth', 1.8);

% System 2: red
plot(real(zeros2), imag(zeros2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.8);
plot(real(poles2), imag(poles2), 'rx', 'MarkerSize', 8, 'LineWidth', 1.8);

title('Pole-Zero Map Comparison');
xlabel('Real Axis');
ylabel('Imaginary Axis');
legend('Unit Circle', 'Zeros (Sys1)', 'Poles (Sys1)', 'Zeros (Sys2)', 'Poles (Sys2)');
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;

% =====================================================
% Plot 2: Frequency Response (0 → π)
% =====================================================
[H1, w] = freqz(b1, a1, 1024, 'half');
[H2, ~] = freqz(b2, a2, 1024, 'half');

figure;
plot(w, abs(H1), 'b', 'LineWidth', 1.8); hold on;
plot(w, abs(H2), 'r', 'LineWidth', 1.8);
grid on;
title('|H(e^{j\omega})| Frequency Response Comparison');
xlabel('Frequency (radians/sample)');
ylabel('Magnitude');
legend('System 1', 'System 2');
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;
