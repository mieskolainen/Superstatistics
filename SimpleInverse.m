% Simple combinatorics forward-inverse closure test
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;
addpath ./src

N = 3;
bitreverse = true; % Bit order makes no difference here
LAMBDA = amat(N, bitreverse);

% Truth
p =[0.040079021973800
    0.050471072739824
    0.222838163465361
    0.173554170276039
    0.218088208991535
    0.164300170447909
    0.130669192105532];

% Poisson mu
R = 0.4;
mu = -log(1 - R)

% Forward map
z = LAMBDA*p;
yhat = inv(LAMBDA) * (exp(-mu*z)-1) / (exp(-mu)-1)

% Inverse map
z = LAMBDA*yhat;
phat = inv(LAMBDA) * log((exp(-mu) - 1)*z + 1) / (-mu)

% Truth - Estimate
dp = p - phat;
stem(dp);
ylabel('$p - \hat{p}$','interpreter','latex');
xlabel('ID','interpreter','latex');



