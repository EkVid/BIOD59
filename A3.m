% Part A 

% Q1 a)

% Parameters
m_3 = 1.2;
good_survival_rate = [0.3, 0.7, 0.85, 0.99];
bad_survival_rate = [0.1, 0.3, 0.75, 0.85];

p_0 = (good_survival_rate(1) + bad_survival_rate(1)) / 2;
p_1 = (good_survival_rate(2) + bad_survival_rate(2)) / 2;
p_2 = (good_survival_rate(3) + bad_survival_rate(3)) / 2;
p_3 = (good_survival_rate(4) + bad_survival_rate(4)) / 2;


% initial population
% new borns, 1-yr olds, 2-yr olds, adults
N0 = [5; 7; 15; 35];

% 40 years
years = 40;


% given matrix in question
postbreeding_matrix = [0, 0, m_3 * p_2, m_3 * p_3;
                       p_0, 0, 0, 0;
                       0, p_1, 0, 0;
                       0, 0, p_2, p_3];


% initialzize population matrix to store future population for each
% generation
population_each_gen = zeros(4, years);

% initialize total population vector
total_pop = zeros(1, years);

% initial total population is the sum of N0
total_pop(1) = sum(N0);

% initialize annual growth rate 
annual_growth_rate = zeros(1, years);

% populate the initial population with 
% initial population for each generation
population_each_gen(:, 1) = N0;

for t = 2 : years
    % perform transition matrix multiplication
    population_each_gen(:, t) = postbreeding_matrix * population_each_gen(:, t - 1);
    % sum over each generation's population to get the total population
    total_pop(t) = sum(population_each_gen(:, t));
    % calculate the annual growth rate
    annual_growth_rate(t) = total_pop(t) / total_pop(t-1);
end


% Plot

y = 1 : years;

figure;
% Subplot 1: Population of each generation over time
subplot(3, 1, 1);
plot(y, population_each_gen(1, :)); % new borns
hold on;
plot(y, population_each_gen(2, :)); % 1-yr olds
plot(y, population_each_gen(3, :)); % 2-yr olds
plot(y, population_each_gen(4, :)); % adults
xlabel('Year');
ylabel('Population Size');
title('Population Size by Generation Over Time');
legend('Newborns (n_0)', '1-yr olds (n_1)', '2-yr olds (n_2)', 'Adults (n_3)');
hold off;

% Subplot 2: Total population size over time
subplot(3, 1, 2);
plot(y, total_pop);
ylabel('Total Population Size');
xlabel('Year');
title('Total Population Size Over Time');

% Subplot 3: Annual growth rate over time
subplot(3, 1, 3);
plot(y, annual_growth_rate);
ylabel('Annual Growth Rate');
xlabel('Year');
title('Annual Growth Rate (N_t / N_{t-1}) Over Time');



% Q1 b)

% Get the eigenvalues and eigenvectors of the postbreeding matrix
[eig_vec, eig_val] = eig(postbreeding_matrix);

% Get the dominant/largest eigenvalue from the matrix returned above
lambda = max(diag(eig_val));

% Get the eigenvector corresponding to dominant eigenvalue, which is 
% stable stage distribution

% First, get the index of the eigen vector
idx = find(diag(eig_val) == lambda);

% Extract the eigenvector
stable_stage_dist = eig_vec(:, idx);

% normalize it
stable_stage_dist = stable_stage_dist / sum(stable_stage_dist);


% Q2 a)

% assume only good years
p_0 = good_survival_rate(1);
p_1 = good_survival_rate(2);
p_2 = good_survival_rate(3);
p_3 = good_survival_rate(4);

% initialize matrices and vectors, same way as above, just override the
% values
postbreeding_matrix = [0, 0, m_3 * p_2, m_3 * p_3;
                       p_0, 0, 0, 0;
                       0, p_1, 0, 0;
                       0, 0, p_2, p_3];

population_each_gen = zeros(4, years);
total_pop = zeros(1, years);
total_pop(1) = sum(N0);
annual_growth_rate = zeros(1, years);
population_each_gen(:, 1) = N0;

for t = 2 : years
    population_each_gen(:, t) = postbreeding_matrix * population_each_gen(:, t - 1);
    total_pop(t) = sum(population_each_gen(:, t));
    annual_growth_rate(t) = total_pop(t) / total_pop(t-1);
end

figure;

% Subplot for Total Population (Good Years)
subplot(2, 1, 1);
plot(1:years, total_pop, 'LineWidth', 2);
xlabel('Year');
ylabel('Total Population Size');
title('Total Population Size (Good Years)');

% Subplot for Annual Growth Rate (Good Years)
subplot(2, 1, 2);
plot(1:years, annual_growth_rate, 'LineWidth', 2);
xlabel('Year');
ylabel('Annual Growth Rate');
title('Annual Growth Rate (Good Years)');

sgtitle('Population Dynamics Under Good Years');


% assume only bad years
p_0 = bad_survival_rate(1);
p_1 = bad_survival_rate(2);
p_2 = bad_survival_rate(3);
p_3 = bad_survival_rate(4);

% initialize matrices and vectors, same way as above, just override the
% values
postbreeding_matrix = [0, 0, m_3 * p_2, m_3 * p_3;
                       p_0, 0, 0, 0;
                       0, p_1, 0, 0;
                       0, 0, p_2, p_3];

population_each_gen = zeros(4, years);
total_pop = zeros(1, years);
total_pop(1) = sum(N0);
annual_growth_rate = zeros(1, years);
population_each_gen(:, 1) = N0;

for t = 2 : years
    population_each_gen(:, t) = postbreeding_matrix * population_each_gen(:, t - 1);
    total_pop(t) = sum(population_each_gen(:, t));
    annual_growth_rate(t) = total_pop(t) / total_pop(t-1);
end


% plot?
figure;
 
% Subplot for Total Population (Bad Years)
subplot(2, 1, 1);
plot(1:years, total_pop, 'LineWidth', 2);
xlabel('Year');
ylabel('Total Population Size');
title('Total Population Size (Bad Years)');

% Subplot for Annual Growth Rate (Good Years)
subplot(2, 1, 2);
plot(1:years, annual_growth_rate, 'LineWidth', 2);
xlabel('Year');
ylabel('Annual Growth Rate');
title('Annual Growth Rate (Bad Years)');

sgtitle('Population Dynamics Under Bad Years');




% Q2 b)

% Parameters

p = 0.5;
simulations = 100;
years = 40;

% matrices to store values
total_pop_sim = zeros(simulations, years);
population_each_gen_sim = zeros(4, years, simulations);

for i = 1 : simulations

    % initialize population matrix
    population_each_gen = zeros(4, years);

    % populate the initial population into the matrix
    population_each_gen(:, 1) = N0;

    for t = 2 : years
        % if random number < 0.5, then I call it good year
        is_good_year = rand < p;

        if is_good_year
            survival_rates = good_survival_rate;
        else
            survival_rates = bad_survival_rate;
        end
        
        % Update postbreeding matrix with current survival rates
        p_0 = survival_rates(1);
        p_1 = survival_rates(2);
        p_2 = survival_rates(3);
        p_3 = survival_rates(4);
        m_3 = 1.2; % fecundity rate

        postbreeding_matrix = [0, 0, m_3 * p_2, m_3 * p_3;
                               p_0, 0, 0, 0;
                               0, p_1, 0, 0;
                               0, 0, p_2, p_3];

        % perform matrix multiplication
        population_each_gen(:, t) = postbreeding_matrix * population_each_gen(:, t - 1);
    end
    % store to the matrics
    population_each_gen_sim(:, :, i) = population_each_gen;
    total_pop_sim(i, :) = sum(population_each_gen, 1);
end

figure;
% Plot the first 9 runs
for i = 1 : 9
    subplot(3, 3, i); % Create a 3x3 grid of subplots
    y = 1 : years;
    
    % Plot total population size
    plot(y, total_pop_sim(i, :), 'DisplayName', 'Total Population Size', 'LineWidth', 1.5);
    hold on;
    
    % Plot generation populations
    plot(y, squeeze(population_each_gen_sim(1, :, i)), 'DisplayName', 'Newborns (n_0)', 'LineWidth', 1);
    plot(y, squeeze(population_each_gen_sim(2, :, i)), 'DisplayName', '1-yr olds (n_1)', 'LineWidth', 1);
    plot(y, squeeze(population_each_gen_sim(3, :, i)), 'DisplayName', '2-yr olds (n_2)', 'LineWidth', 1);
    plot(y, squeeze(population_each_gen_sim(4, :, i)), 'DisplayName', 'Adults (n_3)', 'LineWidth', 1);
    
    % Labels and title
    xlabel('Year');
    ylabel('Population Size');
    title(['Simulation ' num2str(i)]);
    legend('show', 'Location', 'best'); % Show legends for better understanding
    hold off;
end

% find how many runs went to quasi-extinction

% initial population size is 62
N0_total = 62;

% count how many went to quasi extinction
quasi_extinct_count = sum(total_pop_sim(:, end) < 0.2 * N0_total);

% calculate proportion
quasi_extinct_proportion = quasi_extinct_count / simulations;

disp(['Proportion of runs resulting in quasi-extinction: ' num2str(quasi_extinct_proportion)]);


% Bonus

% c)

% redo the same process just update p to 0.25, 0.15, 0.05
% graphs are included in the word doc, avoid redundant copy paste code here

% d) 

% initialize p range
p_range = 0 : 0.01 : 1;

% initialize array/vector for quasi extinction probability
quasi_extinct_prob = zeros(1, length(p_range));

for i = 1 : length(p_range)
    % initialize population matrix
    p = p_range(i);
    total_pop_sim = zeros(simulations, years);

    % run simulations
    for x = 1 : simulations
        population_each_gen = zeros(4, years);
        population_each_gen(:, 1) = N0;

        for t = 2 : years
            % see if it's good or bad year
            is_good_year = rand < p;

            % set values based on good or bad year
            if is_good_year
                survival_rates = good_survival_rate;
            else
                survival_rates = bad_survival_rate;
            end
            
            % Update postbreeding matrix with current survival rates
            p_0 = survival_rates(1);
            p_1 = survival_rates(2);
            p_2 = survival_rates(3);
            p_3 = survival_rates(4);
            m_3 = 1.2; % fecundity rate
    
            postbreeding_matrix = [0, 0, m_3 * p_2, m_3 * p_3;
                                   p_0, 0, 0, 0;
                                   0, p_1, 0, 0;
                                   0, 0, p_2, p_3];

            % perform matrix multiplication
            population_each_gen(:, t) = postbreeding_matrix * population_each_gen(:, t - 1);
        end
         % Store total population size
        total_pop_sim(x, :) = sum(population_each_gen, 1);
    end
    % Quasi-extinction Calculation

    % Count quasi-extinction cases
    quasi_extinction_count = sum(total_pop_sim(:, end) < N0_total * 0.2); 
    
    % Calculate the probability of quasi-extinction
    quasi_extinct_prob(i) = quasi_extinction_count / simulations; 
end

% Plotting the results
figure;
plot(p_range, quasi_extinct_prob, 'LineWidth', 2);
xlabel('Probability of Good Year (p)');
ylabel('Quasi-extinction Probability');
title('Quasi-extinction Probability as a Function of p');