clc, clear, close all
options = optimset('Display', 'off'); x0 = [];
%% Problem 2 Part 2: Repeat Problem 2 Part 1b with More Assets
% Given Information
total_months = (2022-2014)*12; % Number of months elapsed in the data

SPY = readtable('SPY.csv', 'ReadVariableNames', false);   % Import the SPY monthly dataset from Jan 2014 to Jan 2022
GOVT = readtable('GOVT.csv', 'ReadVariableNames', false); % Import the GOVT monthly dataset from Jan 2014 to Jan 2022
EEMV = readtable('EEMV.csv', 'ReadVariableNames', false); % Import the EEMV monthly dataset from Jan 2014 to Jan 2022
CME = readtable('CME.csv', 'ReadVariableNames', false); % Import the CME monthly dataset from Jan 2014 to Jan 2022
BR = readtable('BR.csv', 'ReadVariableNames', false); % Import the BR monthly dataset from Jan 2014 to Jan 2022
CBOE = readtable('CBOE.csv', 'ReadVariableNames', false); % Import the CBOE monthly dataset from Jan 2014 to Jan 2022
ICE = readtable('ICE.csv', 'ReadVariableNames', false); % Import the ICE monthly dataset from Jan 2014 to Jan 2022
ACN = readtable('ACN.csv', 'ReadVariableNames', false); % Import the ACN monthly dataset from Jan 2014 to Jan 2022

% Calculate the Expected Return of each asset
% Find the column of interest and convert table values to numerical format
SPY_adj_close = table2array(SPY(:, 6));
GOVT_adj_close = table2array(GOVT(:, 6));
EEMV_adj_close = table2array(EEMV(:, 6));
CME_adj_close = table2array(CME(:, 6));
BR_adj_close = table2array(BR(:, 6));
CBOE_adj_close = table2array(CBOE(:, 6));
ICE_adj_close = table2array(ICE(:, 6));
ACN_adj_close = table2array(ACN(:, 6));

% Calculate the expected monthly returns of each asset
r_SPY_m = zeros(total_months, 1);
r_GOVT_m = zeros(total_months, 1);
r_EEMV_m = zeros(total_months, 1);
r_CME_m = zeros(total_months, 1);
r_BR_m = zeros(total_months, 1);
r_CBOE_m = zeros(total_months, 1);
r_ICE_m = zeros(total_months, 1);
r_ACN_m = zeros(total_months, 1);
for i = 1:total_months
    r_SPY_m(i) = SPY_adj_close(i+1)/SPY_adj_close(i)-1;
    r_GOVT_m(i) = GOVT_adj_close(i+1)/GOVT_adj_close(i)-1;
    r_EEMV_m(i) = EEMV_adj_close(i+1)/EEMV_adj_close(i)-1;
    r_CME_m(i) = CME_adj_close(i+1)/CME_adj_close(i)-1;
    r_BR_m(i) = BR_adj_close(i+1)/BR_adj_close(i)-1;
    r_CBOE_m(i) = CBOE_adj_close(i+1)/CBOE_adj_close(i)-1;
    r_ICE_m(i) = ICE_adj_close(i+1)/ICE_adj_close(i)-1;
    r_ACN_m(i) = ACN_adj_close(i+1)/ACN_adj_close(i)-1;
end

% Combined monthly returns of each asset
r_m_total = [r_SPY_m, r_GOVT_m, r_EEMV_m, r_CME_m, r_BR_m, r_CBOE_m, r_ICE_m, r_ACN_m];

r_SPY = mean(r_SPY_m);   % Arithmetic average monthly return of $SPY from 2014-2022
r_GOVT = mean(r_GOVT_m); % Arithmetic average monthly return of $GOVT from 2014-2022
r_EEMV = mean(r_EEMV_m); % Arithmetic average monthly return of $EEMV from 2014-2022
r_CME = mean(r_CME_m);   % Arithmetic average monthly return of $CME from 2014-2022
r_BR = mean(r_BR_m);     % Arithmetic average monthly return of $BR from 2014-2022
r_CBOE = mean(r_CBOE_m); % Arithmetic average monthly return of $CBOE from 2014-2022
r_ICE = mean(r_ICE_m);   % Arithmetic average monthly return of $ICE from 2014-2022
r_ACN = mean(r_ACN_m);   % Arithmetic average monthly return of $ACN from 2014-2022
mu = [r_SPY, r_GOVT, r_EEMV, r_CME, r_BR, r_CBOE, r_ICE, r_ACN]; % Vector of expected returns of all assets

% Calculate the covariance of each asset
cov_matrix = zeros(length(mu), length(mu)); % Initialize the total covariance matrix in order of: SPY, GOVT, EEMV, CME, BR, CBOE, ICE, ACN
for i = 1:length(mu)
    for j = 1:length(mu)
        cov_ij = cov(r_m_total(:, i), r_m_total(:, j));
        cov_matrix([i, j], [i, j]) = cov_ij;
    end
end

% Define the goal return of the portfolio
n_points = 10; % Number of equally spaced out points on the efficient frontier
goal_R = min(mu):(max(mu)-min(mu))/(n_points-1):max(mu); % Expected return goals range from minimum asset return to maximum asset return

% Perform quadprog
Q = cov_matrix;
c = zeros(length(mu), 1);
A = -mu;
Aeq = ones(1, length(mu)); beq = 1;
ub = [];
lb = []; % with short selling assumed

std_portfolio = zeros(n_points, 1); % Initialize a vector of standard deviations for each goal R
optimal_weights = zeros(n_points, length(mu)); % Initialize a matrix of weights in each asset for each goal R
for i = 1:n_points
    b = -goal_R(i); % The Markowitz return constraint changes for each goal R
    [x, fval] = quadprog(Q, c, A, b, Aeq, beq, lb, ub, x0, options); % Find weights and 1/2 portfolio variance (objective function)
    std_portfolio(i) = (fval*2)^0.5; % Store the portfolio standard deviation for plotting
    optimal_weights(i, :) = x; % Store the optimal weights for table creation
end

% Plot the efficient frontier
plot(std_portfolio, goal_R, '-ko')
xlabel('Volatility \sigma')
ylabel('Expected Return Goal R')
title('The Efficient Frontier of MVO - Problem 2 Part 2')
grid on

% Print a Table of optimal weights and portfolio variance for each goal R
fprintf('Problem 2 Part 2 Answers:\nTable of Optimal Weights and Portfolio Variance\n')
fprintf('Return goal R   SPY Weight   GOVT Weight   EEMV Weight   ')
fprintf('CME Weight   BR Weight   CBOE Weight   ICE Weight   ACN Weight   ')
fprintf('Portfolio Variance\n')
for i = 1:n_points
    fprintf('%13.3f%13.3f%14.3f%14.3f%13.3f%12.3f%14.3f%13.3f%13.3f%21.6f\n'...
        , goal_R(i), optimal_weights(i,1), optimal_weights(i,2),...
        optimal_weights(i,3), optimal_weights(i,4), optimal_weights(i,5),...
        optimal_weights(i,6), optimal_weights(i,7), optimal_weights(i,8),...
        std_portfolio(i)^2)
end