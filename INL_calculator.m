%% ADC INL estimation via noisy code-averaging and local linear fits
% Author: (you)
% Date: today

clear; clc; close all;
rng(42);                         % reproducibility

%% 1) ADC basics
N  = 10;                         % bits
K  = 2^N;                        % codes: 0..K-1
LSB = 1;                         % since input is in integer LSB units
xmin = 0;
xmax_allowed = K - 0.5;          % 1023.5 for 10-bit

%% 2) Define INL on transitions: 1-cycle sine over full code range, amplitude 0.5 LSB
A_inl = 0.3 * LSB;               % amplitude in LSB
% Transitions are between code n and n+1 at ideal (n+0.5). There are K-1 internal transitions.
n_tr = 0:(K-2);                  % transition indices (for thresholds at n+0.5)
inl_true = A_inl * sin( 2*pi*(n_tr + 0.5)/K );    % one full cycle across code range
T_true = (n_tr + 0.5) + inl_true;                 % true transition levels (LSB units)

%% 3) Create equally-spaced inputs, 3 per LSB (total length = 3*K)
samples_per_LSB = 3;
P = samples_per_LSB * K;         % required length
x_grid = (0:(P-1)) / samples_per_LSB;   % step = 1/3 LSB, spans [0, K - 1/3]
% Keep inputs within allowed ADC range [0, 1023.5]; clamp any above
x_grid = min(x_grid, xmax_allowed);

%% 4) Add Gaussian noise (σ = 0.1 LSB) for M=4 independent repeats
M = 2048;
sigma = 0.1 * LSB;
X = repmat(x_grid, M, 1) + sigma * randn(M, P);   % size M x P

%% 5) Quantize each sample using the true (INL-shifted) transitions
% Use bin edges [-Inf, T_true, +Inf] to discretize -> codes 0..K-1
edges = [-Inf, T_true, Inf];     % 1 x (K+1)
codes = discretize(X, edges) - 1; % M x P, returns 0..K-1 (since discretize bins are 1-based)

%% 6) For each input point, average the M codes
y_avg = mean(double(codes), 1);  % 1 x P

%% 7) Slide L=4-point window, fit y = a*x + b, and find x where y = n + 0.5
Lfit = 4;
cand = cell(1, K-1);             % collect candidate transition x for each transition n=0..K-2

for i = 1:(P - Lfit + 1)
    xi = x_grid(i:i+Lfit-1);
    yi = y_avg(i:i+Lfit-1);

    % Robust linear fit
    p = polyfit(xi, yi, 1);      % yi ≈ p(1)*xi + p(2)
    a = p(1); b = p(2);
    if abs(a) < 1e-12, continue; end

    % The window's y-range from the fitted line (use endpoints to avoid noise wiggle)
    y_end = polyval(p, [xi(1), xi(end)]);
    y_min = min(y_end); y_max = max(y_end);

    % Which transitions (n+0.5) can lie within this window's y-range?
    nmin = ceil(y_min - 0.5);
    nmax = floor(y_max - 0.5);
    nmin = max(nmin, 0);
    nmax = min(nmax, K-2);

    if nmin <= nmax
        for n = nmin:nmax
            x_star = ( (n + 0.5) - b ) / a;      % solve y = n+0.5
            if x_star >= xi(1) && x_star <= xi(end)
                cand{n+1}(end+1) = x_star; %#ok<SAGROW>
            end
        end
    end
end

% Average candidates per transition
T_est = nan(1, K-1);
for n = 0:(K-2)
    if ~isempty(cand{n+1})
        T_est(n+1) = mean(cand{n+1});
    end
end

% If any edge transitions were missed (rare), fill by linear interpolation (no extrapolation)
miss = isnan(T_est);
if any(miss)
    known_idx = find(~miss);
    if numel(known_idx) >= 2
        T_est(miss) = interp1(known_idx, T_est(known_idx), find(miss), 'linear', 'extrap');
    else
        error('Too few transition estimates found — increase sampling/overlap or reduce noise.');
    end
end

%% 8) Calculate measured INL (transition-centric)
inl_est = T_est - (n_tr + 0.5);  % same definition as true INL on transitions

%% 9) Compare measured INL vs injected INL
err = inl_est - inl_true;
max_abs_err = max(abs(err));
rmse = sqrt(mean(err.^2));

fprintf('Max |INL_error| = %.6f LSB\n', max_abs_err);
fprintf('RMSE(INL_error) = %.6f LSB\n', rmse);

%% Plots
figure('Color','w'); hold on; grid on;
plot(n_tr + 0.5, inl_true,  'LineWidth',1.5);
plot(n_tr + 0.5, inl_est,   '--', 'LineWidth',1.5);
xlabel('Ideal transition level (LSB units)');
ylabel('INL (LSB)');
title(sprintf('Injected vs Estimated INL (N=%d, M=%d, L=%d, A=%.1f LSB, \\sigma=%.1f LSB)', ...
      N, M, Lfit, A_inl/LSB, sigma/LSB));
legend('True INL','Estimated INL','Location','best');

% (Optional) show error
figure('Color','w'); plot(n_tr + 0.5, err, 'LineWidth',1.2); grid on;
xlabel('Ideal transition level (LSB units)');
ylabel('INL estimation error (LSB)');
title(sprintf('INL Estimation Error | Max=%.4g LSB, RMSE=%.4g LSB', max_abs_err, rmse));
