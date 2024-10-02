% The Ricker Model
function [pop_size] = Ricker(initial_pop, carrying_capacity, growth_rate)
   next_year_pop = initial_pop;
   % compute the exponent 
   exponent = growth_rate * (1 - next_year_pop/carrying_capacity);
   % compute the final result
   next_year_pop = next_year_pop * exp(exponent);
   pop_size = next_year_pop;
end

% The Beverton and Holt Model
function [pop_size] = B_and_H(initial_pop, carrying_capacity, repro_num)
   numerator = repro_num * initial_pop;
   denom = 1 + (repro_num - 1) * initial_pop / carrying_capacity;
   pop_size = numerator / denom;
end

% function to solve diff eq
function [t, N] = solve_diff_eq(N0, r, K, A, tspan, e)
    % Define the differential equation as an anonymous function for quick
    % calculation
    ode = @(t, N) r * N * (1 - N / K) * (N / A - 1) - e * N;
    
    % Solve the ODE using ode23s
    [t, N] = ode23s(ode, tspan, N0);
end



% Part A Q1
initial_population = 10;
carrying_capacity = 1000;

% Initialize empty matrix with 4 rows and 100 entries each
population_size = zeros(4,100);

% Initialize a custom growth rate array
growth_rates = [1.8, 2.1, 2.6, 3];

% Compute population sizes
for t = 1 : 4
    % Initialize the first year population to be 10 for each growth rate
    population_size(t, 1) = initial_population;
    % Compute population for the next 100 years
    for x = 1 : 100
       population_size(t,x + 1) = Ricker(population_size(t,x), carrying_capacity, growth_rates(t));
    end
end

% Create a new figure for Ricker model plots
figure;
% Generate 4 subplots, one for each growth rate
for t = 1 : 4
    subplot(2, 2, t);  % Create a 2x2 grid of subplots
    plot(0:100, population_size(t, :));  
    title(['Growth Rate = ', num2str(growth_rates(t))]);  
    xlabel('Years');  
    ylabel('Population Size');  
end

% Part A Q2

% Initialize empty matrix with 2 rows and 100 entries each
population_size = zeros(2,100);

% Initialize custom reproductive number R array
reproductive_num = [2, 11];

% Compute population sizes
for t = 1 : 2
    % Initialize the first year population to be 10 for each reproductive number
    population_size(t, 1) = initial_population;
    % Compute population for the next 100 years
    for x = 1 : 100
        population_size(t, x+1) = B_and_H(population_size(t, x), carrying_capacity, reproductive_num(t));
    end
end

% Create a new figure for Beverton and Holt model plots
figure;
% Generate 2 subplots, one for each reproductive number
for t = 1 : 2
    subplot(2, 1, t);  % Create a 2x1 grid of subplots
    plot(0:100, population_size(t, :));  
    title(['Reproductive Number = ', num2str(reproductive_num(t))]); 
    xlabel('Years');  
    ylabel('Population Size');  
end
    
% Bonus 


% B-H

% Define the x range and step
x = 0 : 1 : 200;

% Define the function for B-H
% choose K = 100, r = 2 and simplified
y = (200 * x) ./ (100 + x);

% Diagonal line y2 = x
y2 = x;

% Initial point for cobwebbing
N0 = 10;

% Number of iterations for cobwebbing
iterations = 20;  

% Initialize cobweb points
cobweb_x = N0;
cobweb_y = 0;

% cobwebbing 
current_x = N0;
for i = 1 : iterations
    next_y = (200 * current_x) ./ (100 + current_x);  % Evaluate the function
    cobweb_x = [cobweb_x, current_x, current_x];  % Move horizontally
    cobweb_y = [cobweb_y, current_x, next_y];  % Move vertically
    cobweb_x = [cobweb_x, next_y];  % Move horizontally to y = x line
    cobweb_y = [cobweb_y, next_y];
    current_x = next_y;  
end

% Plot the functions
figure;
plot(x, y, 'b-', 'DisplayName', 'y = (200 * x) ./ (100 + x)'); 
hold on; 
plot(x, y2, 'r-', 'DisplayName', 'y2 = x'); 

% Plot the cobweb
plot(cobweb_x, cobweb_y, 'k-', 'DisplayName', 'Cobweb Path');

% Highlight the initial point
plot(N0, 0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Initial Population');

% labels and legend
title('Cobwebbing on B-H Model');
xlabel('x');
ylabel('y');
legend;
grid on; 
hold off;


% Ricker

% Define the x range and step
x = 0 : 1 : 200;

% Define the function for Ricker
% choose K = 100, r = 2 and simplified
y = x .* exp((100 - x) / 50);

% Diagonal line y2 = x
y2 = x;

% Initial point for cobwebbing
N0 = 10;
iterations = 50;  % Number of iterations for cobwebbing

% Initialize cobweb points
cobweb_x = N0;
cobweb_y = 0;

% cobwebbing
current_x = N0;
for i = 1 : iterations
    next_y = current_x * exp((100 - current_x) / 50);  % Evaluate the function
    cobweb_x = [cobweb_x, current_x, current_x];  % Move horizontally
    cobweb_y = [cobweb_y, current_x, next_y];  % Move vertically
    cobweb_x = [cobweb_x, next_y];  % Move horizontally to y = x line
    cobweb_y = [cobweb_y, next_y];
    current_x = next_y;  
end

% Plot the functions
figure;
plot(x, y, 'b-', 'DisplayName', 'y = x .* exp((100-x)/50)'); 
hold on; 
plot(x, y2, 'r-', 'DisplayName', 'y2 = x'); 

% Plot the cobweb
plot(cobweb_x, cobweb_y, 'k-', 'DisplayName', 'Cobweb Path');

% Highlight the initial point
plot(N0, 0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Initial Population');

title('Cobwebbing on Ricker Model');
xlabel('x');
ylabel('y');
legend;
grid on; 
hold off;



% Part B 

% a)

% Parameters
r = 1;        
K = 1000;     
A = 50;       
tspan = [0, 20]; 

% Initial conditions
N0 = [25, 50, 500, 1500];

figure; 
for i = 1 : length(N0)
   
    % solve the diff eq when e = 0
    [t, N] = solve_diff_eq(N0(i), r, K, A, tspan, 0);
    
    % Create a subplot for each initial condition
    subplot(2, 2, i);  
    plot(t, N, 'DisplayName', ['N0 = ' num2str(N0(i))]);
    xlabel('Time');
    ylabel('Fish Population (tons)');
    title(['Population Dynamics for N0 = ' num2str(N0(i))]);
    grid on;
end

sgtitle('Population Dynamics of Fish Stock with different initial population');


% b)

% 100 points between 0 and 8
e = linspace(0, 8, 10000);

% initialize empty arrays
N_1 = NaN(size(e));
N_2 = NaN(size(e));
N_3 = NaN(size(e));

% simplified equilibria equation when set to 0
% (-r/AK)N^3 + (r/A + r/K)N^2 - (r + e)N = 0


for i = 1 : length(e)
    % (-r/AK)N^3 + (r/A + r/K)N^2 - (r + e)N = 0

    coef = [-r / (A * K), (r / A + r / K), -r - e(i), 0];

    % solve the cubic equation:
    try
        root = roots(coef);
        
        % Only keep real roots
        real_roots = root(imag(root) == 0);
        
        % Sort real roots in ascending order
        real_roots = sort(real_roots);
        
        % Store roots
        if length(real_roots) > 2
            N_1(i) = real_roots(1);
            N_2(i) = real_roots(2);
            N_3(i) = real_roots(3);
        elseif length(real_roots) == 2
            N_1(i) = real_roots(1);
            N_2(i) = real_roots(2);
            N_3(i) = NaN;
        elseif isscalar(real_roots)
            N_1(i) = real_roots(1);
            N_2(i) = NaN;
            N_3(i) = NaN;
        else
            N_1(i) = NaN;
            N_2(i) = NaN;
            N_3(i) = NaN;
        end
    catch
        % If there's an error solving roots, continue to next iteration
        disp(['Error solving for roots with e = ', num2str(e(i))]);
        continue;
    end
end


% plot the equilibria as a function of e
figure;
plot(e, N_1, 'b', 'DisplayName', 'Equilibrium 1');
hold on;
plot(e, N_2, 'r', 'DisplayName', 'Equilibrium 2');
hold on;
plot(e, N_3, 'g', 'DisplayName', 'Equilibrium 3');
xlabel('Harvesting Effort (e)');
ylabel('Equilibrium Population (N)');
title('Equilibrium Population as a Function of Harvesting Effort');
legend show;



% c)

% e for part c
e_c = [3, 4, 4.5125, 5]; 

% N0 for part c
N0_c = [300, 525, 700];  

% time span for part c
tspan_c = [0 5];

figure;

% Loop through each e value
for e_idx = 1 : length(e_c)

    % Loop through each N0
    for N0_idx = 1 : length(N0_c)

        % solve the diff eq based on each N0 value
        [t, N] = solve_diff_eq(N0_c(N0_idx), r, K, A, tspan_c, e_c(e_idx));

        % Create a subplot for each case
        subplot(length(e_c), length(N0_c), (e_idx - 1) * length(N0_c) + N0_idx);
        plot(t, N, 'LineWidth', 2);
        xlabel('Time');
        ylabel('Fish Population (tons)');
        ylim([0, 800]);
        title(['e = ' num2str(e_c(e_idx)) ', N0 = ' num2str(N0_c(N0_idx))]);
    end
end

sgtitle('Population Dynamics of Fish Stock with Different Harvesting Efforts and Initial Populations');


