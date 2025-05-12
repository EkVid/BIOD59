
% load the data downloaded from Quercus
p = struct();
run("BIOD59_HW4_Solution_Parameters.m");

% a)

% idea: use polyfit to find a linear equation like y = mx + b to fit all
% the log transformted data (y) and (1/T - 1/T0) (x)
% then use the returned [m,b] and polyval to calculate the exact values at (1/T - 1/T0)
% plot the points ( x, y ) and plot the linear regression fit by using the
% points obtained from polyval

% initialize vectors for storage of activation energies and data's
% reference rate for theta, mu, and rho
Es = zeros(1,3);
D_0s = zeros(1, 3);

% names and greek letters of the processes
process_names = {'Development', 'Mortality', 'Transmission'};
symbols = ["θ_0", "μ_0", "ρ_0"];
E_symbols = ["E_θ", "E_μ", "E_ρ"];

% figure for plotting purpose
figure;

% loop over each thermal performance dataset
for i = 1 : 3
    % Log Transform the data
    log_D = log(p.TP_Data(i, :));

    % Linear Regression
    % coef(1) is the slope
    % coef(2) is the intercept for y = mx + b
    coef = polyfit((1 ./p.T_Data - 1 / p.T0), log_D, 1);

    % stores E_theta, E_mu, and E_rho into the vector
    Es(i) = -coef(1) * p.k;

    % store the original rate parameter (theta_0, mu_0, and rho_0) before log transformation.
    % We need to exponentiate the intercepts in order to get that
    D_0s(i) = exp(coef(2));

    % Print each value out for the 3 equations
    fprintf('%s: %s = %.4f, %s = %.4f\n', process_names{i}, symbols{i}, D_0s(i), E_symbols{i}, Es(i));

    % evaluate them at the points
    pred_log_D = polyval(coef, 1 ./p.T_Data - 1 / p.T0);

    % plot the results
    subplot(3, 1, i);
    % plot the points with x being 1/T - 1/T0, and y being log transformed data 
    plot(1 ./ p.T_Data - 1 / p.T0, log_D, 'bo', 'MarkerSize', 8, 'DisplayName', 'Data');
    hold on;
    % plot the best fit linear regression line
    plot(1 ./ p.T_Data - 1 / p.T0, pred_log_D, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
    xlabel('1/TK - 1/T_{0}K');
    ylabel('log(D(T))');
    title(sprintf('Log-transformed %s Data and Linear Fit', process_names{i}));
    legend('Location', 'best');
    hold off;
end


% b) 


% Initialize vectors for thermal performance at different temperatures 
% using the filled in D_0s and Es
theta_T = D_0s(1) * exp(-Es(1) * (1 ./ (p.T_Vec * p.k) - 1 ./ (p.T0 * p.k)));
mu_T    = D_0s(2) * exp(-Es(2) * (1 ./ (p.T_Vec * p.k) - 1 ./ (p.T0 * p.k)));
rho_T   = D_0s(3) * exp(-Es(3) * (1 ./ (p.T_Vec * p.k) - 1 ./ (p.T0 * p.k)));


% calculate R_0(T) based on the equation given
R0_T = (theta_T ./ (mu_T + theta_T)) .* (rho_T * p.H ./ (mu_T + rho_T * p.H)) .* (p.lambda / (p.mu_P + p.b + p.alpha));


% Create figures with four subplots
figure;

subplot(4, 1, 1); % Development Rate theta(T)
% plot the data points with original temperature as x coordinate 
plot(p.T_Data - 273.15, p.TP_Data(1, :), 'bo', 'MarkerSize', 8, 'DisplayName', 'Data');
hold on;
% plot the thermal performance value for the development rate from the
% calculation with the original temperature
plot(p.T_Vec - 273.15, theta_T, 'r-', 'LineWidth', 2, 'DisplayName', 'Best-Fit Curve');
xlabel('Temperature (°C)'); % x label
ylabel('\theta(T)'); % y label
title('Development Rate \theta(T)');
legend('Location', 'best');
hold off;

subplot(4, 1, 2); % Mortality Rate mu(T)
% plot the data points with original temperature as x coordinate 
plot(p.T_Data - 273.15, p.TP_Data(2, :), 'bo', 'MarkerSize', 8, 'DisplayName', 'Data');
hold on;
% same idea as plotting the development rate
plot(p.T_Vec - 273.15, mu_T, 'r-', 'LineWidth', 2, 'DisplayName', 'Best-Fit Curve');
xlabel('Temperature (°C)'); % x label
ylabel('\mu(T)'); % y label
title('Mortality Rate \mu(T)');
legend('Location', 'best');
hold off;

subplot(4, 1, 3); % Transmission Rate rho(T)
% plot the data points with original temperature as x coordinate 
plot(p.T_Data - 273.15, p.TP_Data(3, :), 'bo', 'MarkerSize', 8, 'DisplayName', 'Data');
hold on;
% same idea as plotting the development rate
plot(p.T_Vec - 273.15, rho_T, 'r-', 'LineWidth', 2, 'DisplayName', 'Best-Fit Curve');
xlabel('Temperature (°C)'); % x label
ylabel('\rho(T)'); % y label
title('Transmission Rate \rho(T)');
legend('Location', 'best');
hold off;

subplot(4, 1, 4); % Basic Reproductive Number R0(T)
% plot the basic reproductive number as a function of temperature, 
% calculated based on the thermal performance curves
plot(p.T_Vec - 273.15, R0_T, 'k-', 'LineWidth', 2, 'DisplayName', 'R_0(T)');
xlabel('Temperature (°C)'); % x label
ylabel('R_0(T)'); % y label
title('Basic Reproductive Number R_0(T)');
legend('Location', 'best');

sgtitle('Thermal Performance and Basic Reproductive Number');


% c)

% Host-parasite model based on given differential equations
host_parasite_model = @(t, N, T) [
    % N is a vector containing F, L, and P, N(1) = F, N(2) = L, N(3) = P
    % T is a vector of temperature in Kelvin

    % dF/dt calculation based on the equation given
    p.lambda * N(3) ...
    - (D_0s(2) * exp(-Es(2) * (1 / (T * p.k) - 1 / (p.T0 * p.k)))) * N(1) ...
    - (D_0s(1) * exp(-Es(1) * (1 / (T * p.k) - 1 / (p.T0 * p.k)))) * N(1); 

    % dL/dt calculation based on the equation given
    (D_0s(1) * exp(-Es(1) * (1 / (T * p.k) - 1 / (p.T0 * p.k)))) * N(1) ...
    - (D_0s(2) * exp(-Es(2) * (1 / (T * p.k) - 1 / (p.T0 * p.k)))) * N(2) ...
    - (D_0s(3) * exp(-Es(3) * (1 / (T * p.k) - 1 / (p.T0 * p.k)))) * N(2) * p.H;

    % dP/dt calculation based on the equation given
    (D_0s(3) * exp(-Es(3) * (1 / (T * p.k) - 1 / (p.T0 * p.k)))) * N(2) * p.H ...
    - (p.mu_P + p.b) * N(3) ...
    - p.alpha * p.H * (N(3) / p.H + N(3)^2 * (p.j + 1) / (p.H^2 * p.j)) 
];

% initialize temperatures as requested in the question in Celsius
temp = [5, 15, 25, 35];

% transform the Celsius degree to Kelvins and store it in temp_k
temp_K = temp + 273.15;

% initialize cell to store results
res = cell(1, length(temp_K));

% Loop from 1 to length(temp_K)
for i = 1 : length(temp_K)
    % Solve the ODE for each temperature
    [t, N] = ode23(@(t,N) host_parasite_model(t, N, temp_K(i)), p.tSpan, p.N0);
    % Store time and adult parasite (N(3)) abundance for each temperature
    res{i} = [t, N(:, 3)];
end

% Plotting the results
figure;
hold on;
colors = ['b', 'g', 'r', 'k'];
for i = 1:length(temp_K)
    plot(res{i}(:, 1), res{i}(:, 2), 'Color', colors(i), 'LineWidth', 2, 'DisplayName', sprintf('T = %d°C', temp(i)));
end
xlabel('Time (days)');
ylabel('Adult Parasite Abundance');
title('Host-Parasite Model: Adult Parasite Abundance Over Time');
legend('Location', 'best');
hold off;