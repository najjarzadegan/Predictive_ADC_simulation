% Generate some sample data
x = 2:4096;
y = Xr(2:end);

% Define the model function
model = @(a, x) a(1)*sin(a(2)*x + a(3)) + a(4);

% Fit the model to the data
f = fit(x', y', model, 'StartPoint', [1, 1, 0, 0]);

% Plot the data and the fitted curve
plot(x, y, 'o', x, f(x))