clc, clear, close all
options = optimset('Display', 'off'); x0 = [];
%% Problem 2 Part 1a): Find Expected Return, Standard Deviation, and Covariances of $SPY, $GOVT, and $EEMV
% Given Information
total_months = (2022-2014)*12; % Number of months elapsed in the data

SPY = readtable('SPY.csv', 'ReadVariableNames', false);   % Import the SPY monthly dataset from Jan 2014 to Jan 2022
GOVT = readtable('GOVT.csv', 'ReadVariableNames', false); % Import the GOVT monthly dataset from Jan 2014 to Jan 2022
EEMV = readtable('EEMV.csv', 'ReadVariableNames', false); % Import the EEMV monthly dataset from Jan 2014 to Jan 2022
% Calculate the Expected Return of these Three ETFs
% Find the column of interest and convert table values to numerical format
SPY_adj_close = table2array(SPY(:, 6));
GOVT_adj_close = table2array(GOVT(:, 6));
EEMV_adj_close = table2array(EEMV(:, 6));

% Calculate the expected yearly returns of each ETF
r_SPY_m = zeros(total_months, 1);
r_GOVT_m = zeros(total_months, 1);
r_EEMV_m = zeros(total_months, 1);
for i = 1:total_months
    r_SPY_m(i) = SPY_adj_close(i+1)/SPY_adj_close(i)-1;
    r_GOVT_m(i) = GOVT_adj_close(i+1)/GOVT_adj_close(i)-1;
    r_EEMV_m(i) = EEMV_adj_close(i+1)/EEMV_adj_close(i)-1;
end

r_SPY = mean(r_SPY_m);   % Arithmetic average monthly return of $SPY from 2014-2022
r_GOVT = mean(r_GOVT_m); % Arithmetic average monthly return of $GOVT from 2014-2022
r_EEMV = mean(r_EEMV_m); % Arithmetic average monthly return of $EEMV from 2014-2022

% Calculate the standard deviation of each ETF
sigma_SPY = std(r_SPY_m);
sigma_GOVT = std(r_GOVT_m);
sigma_EEMV = std(r_EEMV_m);

% Calculate the covariance of each ETF
cov_SPY_GOVT_matrix = cov(r_SPY_m, r_GOVT_m);   % 2x2 Matrix [var(spy), cov(spy, govt); cov(govt, spy), var(govt)]
cov_SPY_EEMV_matrix = cov(r_SPY_m, r_EEMV_m);   % 2x2 Matrix [var(spy), cov(spy, eemv); cov(eemv, spy), var(eemv)]
cov_GOVT_EEMV_matrix = cov(r_GOVT_m, r_EEMV_m); % 2x2 Matrix [var(govt), cov(govt, eemv); cov(eemv, govt), var(eemv)]

% Extract each unique individual covariance
var_SPY = cov_SPY_GOVT_matrix(1);
var_GOVT = cov_GOVT_EEMV_matrix(1);
var_EEMV = cov_GOVT_EEMV_matrix(4);

cov_SPY_GOVT = cov_SPY_GOVT_matrix(2);
cov_SPY_EEMV = cov_SPY_EEMV_matrix(2);
cov_GOVT_EEMV = cov_GOVT_EEMV_matrix(2);

% Print Answers for Part a)
fprintf('Problem 2 Part 1a) Answers (Monthly Values):\n')
fprintf('For SPY:\n\t Expected Return: %.8f\n', r_SPY)
fprintf('\t Standard Deviation: %.8f\n', sigma_SPY)
fprintf('\t Variance: %.8f\n', var_SPY)
fprintf('\t Covariance with GOVT: %.8f\n', cov_SPY_GOVT)
fprintf('\t Covariance with EEMV: %.8f\n\n', cov_SPY_EEMV)

fprintf('For GOVT:\n\t Expected Return: %.8f\n', r_GOVT)
fprintf('\t Standard Deviation: %.8f\n', sigma_GOVT)
fprintf('\t Variance: %.8f\n', var_GOVT)
fprintf('\t Covariance with SPY: %.8f\n', cov_SPY_GOVT)
fprintf('\t Covariance with EEMV: %.8f\n\n', cov_GOVT_EEMV)

fprintf('For EEMV:\n\t Expected Return: %.8f\n', r_EEMV)
fprintf('\t Standard Deviation: %.8f\n', sigma_EEMV)
fprintf('\t Variance: %.8f\n', var_EEMV)
fprintf('\t Covariance with SPY: %.8f\n', cov_SPY_EEMV)
fprintf('\t Covariance with GOVT: %.8f\n\n', cov_GOVT_EEMV)

%% Problem 2 Part 1b): Generate an Efficient Frontier of the Three Assets
% Define the goal return of the portfolio
mu = [r_SPY, r_GOVT, r_EEMV]; % Vector of expected returns of all three ETFs
n_points = 10; % Number of equally spaced out points on the efficient frontier
goal_R = min(mu):(max(mu)-min(mu))/(n_points-1):max(mu); % Expected return goals range from minimum asset return to maximum asset return

% Perform quadprog
Q = [sigma_SPY^2, cov_SPY_GOVT, cov_SPY_EEMV;
    cov_SPY_GOVT, sigma_GOVT^2, cov_GOVT_EEMV;
    cov_SPY_EEMV, cov_GOVT_EEMV, sigma_EEMV^2];
c = [0, 0, 0]';
A = -mu;
Aeq = [1, 1, 1]; beq = 1;
ub = [];
lb = []; % with short selling assumed

std_portfolio = zeros(n_points, 1); % Initialize a vector of standard deviations for each goal R
optimal_weights = zeros(n_points, length(mu)); % Initialize a matrix of weights in each ETF for each goal R
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
title('The Efficient Frontier of MVO - Problem 2 Part 1b)')
grid on

% Print a Table of optimal weights and portfolio variance for each goal R
fprintf('Problem 2 Part 1b) Answers:\nTable of Optimal Weights and Portfolio Variance\n')
fprintf('Return goal R   SPY Weight   GOVT Weight   EEMV Weight   Portfolio Variance\n')
for i = 1:n_points
    fprintf('%13.3f%13.3f%14.3f%14.3f%21.6f\n', goal_R(i), ...
        optimal_weights(i,1), optimal_weights(i,2), optimal_weights(i,3), std_portfolio(i)^2)
end

%% Problem 2 Part 1c): Compare Portfolios and Explain Relative Risk and Return
% Define the portfolio weights
min_var_port = optimal_weights(1, :); % Minimum variance optimal portfolio
equal_port = [1/3, 1/3, 1/3]; % Equal weighted portfolio
given_port = [0.7, 0.2, 0.1]; % Portfolio defined in the question (70% SPY, 20% GOVT, 10% EEMV)

% Calculate the Feb 2022 return for each of the Three ETFs
EoFeb_SPY = 435.28; EoJan_SPY = 448.52; % Adjusted end of month closing prices for SPY (Yahoo Finance USD)
EoFeb_GOVT = 25.69; EoJan_GOVT = 25.91; % Adjusted end of month closing prices for GOVT (Yahoo Finance USD)
EoFeb_EEMV = 62.38; EoJan_EEMV = 62.45; % Adjusted end of month closing prices for EEMV (Yahoo Finance USD)

r_Feb_SPY = EoFeb_SPY/EoJan_SPY-1;    % Feb 2022 return for SPY
r_Feb_GOVT = EoFeb_GOVT/EoJan_GOVT-1; % Feb 2022 return for GOVT
r_Feb_EEMV = EoFeb_EEMV/EoJan_EEMV-1; % Feb 2022 return for EEMV
r_Feb_ETF = [r_Feb_SPY; r_Feb_GOVT; r_Feb_EEMV]; % Feb 2022 ETF returns in vector notation

% Calculate the Feb 2022 return for each portfolio
r_Feb_minvar = min_var_port*r_Feb_ETF;
r_Feb_equal = equal_port*r_Feb_ETF;
r_Feb_given = given_port*r_Feb_ETF;

% Print answers for Part c)
fprintf('\nProblem 2 Part 1c) Answers:\nRank the 3 Portfolios by Return\n')
fprintf('\tRank 1: Minimum Variance Portfolio. Feb 2022 Return = %.2f%%\n', r_Feb_minvar*100)
fprintf('\tRank 2: Equal Weighted Portfolio. Feb 2022 Return = %.2f%%\n', r_Feb_equal*100)
fprintf('\tRank 3: 70%% SPY, 20%% GOVT, 10%% EEMV Portfolio. Feb 2022 Return = %.2f%%\n', r_Feb_given*100)
fprintf('\nExplained relative performance in terms of risk and return in written report\n\n')