% =====================================================
% Compare two transfer functions (H1 and H2)
% =====================================================

% ---- Parameters: First system ----
z0_1 = 0.9;
p0_1 = 0.5;
z1_1 = 0.9;
p1_1 = 0.5;
theta_z1 = 0*pi/180;
theta_p1 = 0*pi/180;

% ---- Compute zeros and poles ----
zeros1 = [z0_1, z1_1*exp(-1j*theta_z1), z1_1*exp(1j*theta_z1)];
poles1 = [p0_1, p1_1*exp(-1j*theta_p1), p1_1*exp(1j*theta_p1)];

% ---- Compute coefficients ----
b1 = poly(zeros1);
a1 = poly(poles1);


% =====================================================
% Plot 1: Pole–Zero Map (custom color control)
% =====================================================
figure; hold on; grid on; axis equal;

% Plot unit circle
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);

% ---- Options for colored arcs on unit circle ----
show_arcs = 1;        % set to 1 to draw arcs, 0 to skip
X_deg     = 10;       % end of green arc (0 -> X); red arc is (X -> 180)

if show_arcs
    % Clamp X to [0, 180]
    X_deg = max(0, min(180, X_deg));

    % Build angles (in radians)
    t_green = deg2rad(linspace(0,   X_deg, 200));   % 0 to X
    t_red   = deg2rad(linspace(X_deg, 180, 200));   % X to 180

    % Arc coordinates (unit circle)
    xg = cos(t_green); yg = sin(t_green);
    xr = cos(t_red);   yr = sin(t_red);

    % Draw arcs
    plot(xg, yg, 'g-', 'LineWidth', 3);   % BOI arc (green)
    plot(xr, yr, 'r-', 'LineWidth', 3);   % OOB arc (red)

    % Annotations placed a bit outside the circle
    mid_green = deg2rad((0 + X_deg)/2);
    mid_red   = deg2rad((X_deg + 180)/2);

    r_txt = 1.2;  % radius for text placement
    text(r_txt*cos(mid_green), r_txt*sin(mid_green), 'BOI', ...
        'Color', 'g', 'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');

    text(0.8*cos(mid_red),   0.8*sin(mid_red),   'OOB', ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 16,'HorizontalAlignment', 'center');
end

% System 1: blue
plot(real(zeros1), imag(zeros1), 'bo', 'MarkerSize', 8, 'LineWidth', 1.8);
plot(real(poles1), imag(poles1), 'bx', 'MarkerSize', 8, 'LineWidth', 1.8);


title('Type 2');
xlabel('Real Axis');
ylabel('Imaginary Axis');
%legend('Unit Circle', 'Zeros (Sys1)', 'Poles (Sys1)');
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;

