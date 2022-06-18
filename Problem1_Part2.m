clc, clear, close all
%% Problem 1 Part 2: Formulate a linear program with Maximum 50% Invested in B-rated Bonds
%% Given information
date = [1, 2, 3, 4, 5, 6]; % Given due dates of liabilities
req_amt = [500, 200, 800, 400, 700, 900]; % Amounts owed at each date
spot = [0.01, 0.015, 0.02, 0.025, 0.03, 0.035]; % Given spot rates

% Bond prices in chronological order
bond_prices = [108, 94, 99, 92.7, 96.6, 95.9, 92.9, 110, 104, 101,...
    107, 102, 95.2]; % Unit bond prices for each bond
reinvestment_cost = [0, 0, 0, 0, 0]; % No extra cost to reinvesting remainders at each time period

% Find short forward rates
f = zeros(1, length(date)); % Initialize the short forward rate vector

% Populate short forward rate vector
for i = 1:length(date)
    if i == 1
        f(i) = spot(i);
    else
        f(i) = (1+spot(i))^i/(1+spot(i-1))^(i-1)-1;
    end
end

% Unit bond payments at each date
coupons = [
    10, 7, 8, 6, 7, 6, 5, 10, 8, 6, 10, 7, 100;    % Each unit bond payment at date 1
    10, 7, 8, 6, 7, 6, 5, 10, 8, 6, 110, 107, 0;   % Each unit bond payment at date 2
    10, 7, 8, 6, 7, 6, 5, 110, 108, 106, 0, 0, 0;  % Each unit bond payment at date 3
    10, 7, 8, 6, 7, 106, 105, 0, 0, 0, 0, 0, 0;    % Each unit bond payment at date 4
    10, 7, 8, 106, 107, 0, 0, 0, 0, 0, 0, 0, 0;    % Each unit bond payment at date 5
    110, 107, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0    % Each unit bond payment at date 6
    ];

reinvestment = [
    -1, 0, 0, 0, 0;      % Remainder at date 1 to be reinvested
    1+f(2), -1, 0, 0, 0; % Remainder at date 1 reinvested at f12, remainder at date 2 to be reinvested
    0, 1+f(3), -1, 0, 0; % Remainder at date 2 reinvested at f23, remainder at date 3 to be reinvested
    0, 0, 1+f(4), -1, 0; % Remainder at date 3 reinvested at f34, remainder at date 4 to be reinvested
    0, 0, 0, 1+f(5), -1; % Remainder at date 4 reinvested at f45, remainder at date 5 to be reinvested
    0, 0, 0, 0, 1+f(6)   % Remainder at date 5 reinvested at f56, any remainder at the end is profit
    ];

B_rated_limits = [
    108, 94, 99, 92.7, 96.6, 95.9,...              % Positive on all B-rated bonds
    -92.9, -110, -104, -101, -107, -102, -95.2,... % Negative on all A-rated bonds
    0, 0, 0, 0, 0                                  % No cost for reinvesting
    ];
B_rated_inequality = 0;                            % B rated bond value - A rated bond value <= 0

%% Define the Objective Coefficients and Perform Linprog
c = [bond_prices, reinvestment_cost]'; % Minimize cost of bond portfolio

% Inequality constraints (Negative to adjust for inequality for linprog)
A = [
    -coupons, -reinvestment; % Coupon payments with remainders after obligations reinvested
    B_rated_limits           % Maximum 50 percent B-rated bonds
    ];
b = [-req_amt, B_rated_inequality]'; % Each obligation must be at least met

% Equality constraints
Aeq = [];
beq = [];

% Variable bounds
ub = []';
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';

%% Print linprog results
[x, fval] = linprog(c, A, b, Aeq, beq, lb, ub);
fprintf('Problem 1 Part 2\n')
fprintf('Investing in:\n')
for i = 1:size(coupons, 2)
    fprintf('\t %.4f units of bond %d\n', x(i), i)
end
fprintf('\nAmount reinvested at each time period\n')
for i = 1:size(reinvestment, 2)
    fprintf('\t Date %d: $%.2f\n', i, x(i+size(coupons, 2)))
end
fprintf('\nTotal cost of the bond portfolio: $%.2f\n\n', fval)

B_value = 0;
for i = 1:6
    B_value = B_value + x(i)*bond_prices(i);
end
fprintf('Weight of the portfolio in B-rated bonds: %.2f%%\n\n', B_value/fval*100)