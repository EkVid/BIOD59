% Q1 b)
initial_pop = 1000;
growth_rate = 1.1;
years = 15;

% Initialize a zero array with the size of years
population = zeros(1, years);

% Assign initial population to the first element in the array
population(1) = initial_pop;

% Loop from 1 to 14 to calculate the new population for each year
for t = 1 : years-1
     population(t+1) = growth_rate * population(t);
end

figure;
plot(0:years-1, population);
title("Population Trajectory for the next 15 years");
xlabel("Years");
ylabel("Population");
grid on;



% Q1 c)
initial_population = 1000;
total_num_growth_rates = 12;

% Define 12 different growth rates
custom_growth_rates = [0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.02];

% Initialize an array to store the last year's population
last_year_populations = zeros(1, total_num_growth_rates);

% Calculate last year's population for each growth rate
for x = 1 : total_num_growth_rates
     last_year_populations(x) = computeLastYearPop(initial_population, custom_growth_rates(x));
end

figure;
plot(custom_growth_rates, last_year_populations);
title("Population Size under different Growth Rates at year 15");
xlabel("Lambda");
ylabel("Population");
grid on;




% Q2 b)
initial_pop = 50;
p_success = 0.8;
years = 10;

% Initialize a zero array with the size of years
survivors = zeros(1, years+1);

% Initial population size is 50 people
survivors(1) = initial_pop;

for t = 1 : years
    survivors(t+1) = binornd(survivors(t), p_success);
end

figure;
plot(0:years, survivors);
title("Population Trajectory for survivors for the next 10 years");
xlabel("Years");
ylabel("Population");
grid on;


% function to compute the last year (15th year) population from the defined 
% growth_rate and initial population
function [last_year_pop] = computeLastYearPop(initial_pop, growth_rate)
   next_year_pop = initial_pop;
   for t = 1 : 15
       next_year_pop = next_year_pop * growth_rate;
   end
   last_year_pop = next_year_pop;
end
