% Statistical uncertainty/error estimators via MC
% bootstrapping and comparison with analytical formulas
% for the combinatorial measurement.
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;


%% DATA

% Number of events after selections (crucial)
N234050 = 3.3e6;
N234039 = 6.1e6;

run234050 = [
0.0
0.019225398
0.007212505
0.018189025
0.001567182
0.000343553
0.002978961
0.009264512
0.003430722
0.000393448
0.000290052
0.00095251
0.002168617
0.000668771
0.00514849
0.028911608
0.039237781
0.005341455
0.001964228
0.007030668
0.000754133
0.000309889
0.002068831
0.009579804
0.015055298
0.002944698
0.002232942
0.007505564
0.014727096
0.006712362
0.074700069
0.709089827];

% Data probabilities
run234039 = [
0
0.019798
0.008434
0.018332
0.001846
0.000346
0.003934
0.009326
0.004642
0.000372
0.000295
0.000963
0.002849
0.000646
0.005106
0.028603
0.040032
0.005213
0.001951
0.006999
0.000762
0.000298
0.002060
0.009438
0.015456
0.002956
0.002224
0.007625
0.015613
0.006501
0.074415
0.702963];

%% Take the two data samples for a comparison

% Drop the first zero ID, the rest sum to one
p_a = run234039(2:end);
p_b = run234050(2:end);

% Sample sizes (number of events after selections)
Na = N234039;
Nb = N234050;

% Dimension, number of observables (detectors)
N = log2(length(p_a) + 1);


%% FULLY SIMULATED DATA SET (JUST FOR TESTS)

%{
% Create "golden ground truth probabilities" from uniform distribution
% (no physical meaning). These have no sample size or uncertainty, because
they are ground truth.
p = rand(2^N-1,1); p = p/sum(p);

% Create two instances, which represent our two measurements
S_a = randsample(1:2^N-1, Na, true, p)';
S_b = randsample(1:2^N-1, Nb, true, p)';

% Put into a vector
p_a = zeros(size(p));
p_b = zeros(size(p));
for i = 1:length(p)
    p_a(i) = sum(S_a == i);
    p_b(i) = sum(S_b == i);
end
% Normalize
p_a = p_a / sum(p_a);
p_b = p_b / sum(p_b);
%}

%% MC Bootstrapping routine

% Number of Monte Carlo (Bootstrap) samples (should be at least 1000 to 10000)
N_MC = 1e3;

% Make MC sample matrices
Xa = zeros(2^N-1, N_MC);
Xb = zeros(2^N-1, N_MC);
XR = zeros(2^N-1, N_MC);

for k = 1:N_MC
    
    % Bootstrap MC sample based on the measurement
    S_a = randsample(1:2^N-1, Na, true, p_a)';
    S_b = randsample(1:2^N-1, Nb, true, p_b)';
    
    % Put events into a vector and classify the IDs
    x_a = zeros(size(p_a));
    x_b = zeros(size(p_b));
    for i = 1:length(p_a)
        x_a(i) = sum(S_a == i);
        x_b(i) = sum(S_b == i);
    end
    
    % Normalize to a probability density
    Xa(:,k) = x_a / sum(x_a);
    Xb(:,k) = x_b / sum(x_b);
    
    % Calculate pointwise ratio between probability distributions A/B
    XR(:,k) = Xa(:,k)./Xb(:,k);
    
    if (mod(k,20) == 0)
       fprintf('MC sample %d/%d \n', k, N_MC); 
    end
end

% Calculate MC estimator 1 sigma confidence intervals
% 1 sigma = 68% Gaussian integral -> 0 ... [16 50 84] ... 1 
% First calculate upper and lower 1 sigma percentiles, and take the
% mean => symmetric upper and lower error bars
sigma_a = 0.5*(abs(prctile(Xa', 16)-prctile(Xa',50)) + abs(prctile(Xa', 84)-prctile(Xa',50)) );
sigma_b = 0.5*(abs(prctile(Xb', 16)-prctile(Xb',50)) + abs(prctile(Xb', 84)-prctile(Xb',50)) );
sigma_R = 0.5*(abs(prctile(XR', 16)-prctile(XR',50)) + abs(prctile(XR', 84)-prctile(XR',50)) );


%% Analytical error formulas are based on Binomial proportion approximations
% The simplified estimator
% https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
% and the ratio is based on first order Taylor expansion of the ratio
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas

% Binomial error (1 sigma, i.e. z_(1-alpha/2), alpha = 68% => z = 1.0)
% 2 sigma, i.e. z = 1.96, would be times 1.96 x times below and so on
sigmaBinom = @(p,N) sqrt(p.*(1-p) / N);

% Ratio error f = A/B, with covariance of A and B (rho_ab = cov_ab/sigma_a*sigma_b)
% which can be set to 0 if unknown
sigmaRatio = @(A,B,sigma_a,sigma_b,cov_ab) ...
    sqrt((A./B).^2.*((sigma_a./A).^2 + (sigma_b./B).^2) - 2*(cov_ab)./(A.*B));


%% Plots

delta = 0.2; % Visual shift on x-axis

% Sample A
figure;
errorbar(1:2^N-1, p_a, sigma_a, 'k.'); hold on;            % MC (black)
errorbar([1:2^N-1]+delta, p_a, sigmaBinom(p_a, Na), 'r.'); % Analytic (red)

axis square; axis tight;
xlabel('$i$ (ID)','interpreter','latex');
ylabel('$P_i$','interpreter','latex');
title(sprintf('Run $A$ (234039), sample size $N_A$ = %0.1E', Na),'interpreter','latex');
set(gca,'yscale','log');
set(gca,'xtick', 1:2:2^N-1);

% Sample B
figure;
errorbar(1:2^N-1, p_b, sigma_b, 'k.'); hold on;            % MC (black)
errorbar([1:2^N-1]+delta, p_b, sigmaBinom(p_b, Nb), 'r.'); % Analytic (red)

axis square; axis tight;
xlabel('$i$ (ID)','interpreter','latex');
ylabel('$P_i$','interpreter','latex');
title(sprintf('Run $B$ (234050), sample size $N_B$ = %0.1E', Nb),'interpreter','latex');
set(gca,'yscale','log');
set(gca,'xtick', 1:2:2^N-1);

% Ratio A/B
figure;
plot(1:2^N-1, ones(2^N-1,1), 'k--'); hold on;
errorbar(1:2^N-1, p_a ./ p_b , sigma_R, 'k.');                            % MC (black)
errorbar([1:2^N-1]+delta, p_a ./ p_b, ...
    sigmaRatio(p_a,p_b,sigmaBinom(p_a, Na),sigmaBinom(p_b, Nb), 1*sigmaBinom(p_a, Na).*sigmaBinom(p_b, Nb)), 'r.'); % Analytic (red)

axis square; axis tight;
title('Black: Bootstrap (MC) $\pm1\sigma$ uncert., Red: Analytic $\pm1\sigma$ uncert.','interpreter','latex');
xlabel('$i$ (ID)','interpreter','latex');
ylabel('$P_i^A/P_i^B \;\;\; (234039/234050)$','interpreter','latex');
set(gca,'xtick', 1:2:2^N-1);

% (A-B)/sigma(B) ~ significance
figure;
plot(1:2^N-1, zeros(2^N-1,1), 'k--'); hold on;
stem(1:2^N-1, (p_a-p_b)./sigmaBinom(p_a, Nb), 'k.', 'linewidth', 2.5);
%stem(1:2^N-1, (p_a*Nb-p_b*Nb)./sqrt(p_b*Nb), 'b.', 'linewidth', 2.5); %
%Poisson approximation

axis square; axis tight;
xlabel('$i$ (ID)','interpreter','latex');
ylabel('$\sim$Significance $(P_i^A-P_i^B)/\sigma_i^B \;\;\; (234039/234050)$','interpreter','latex');
set(gca,'xtick', 1:2:2^N-1);


%% Correlation matrices

% Correlation coefficient of multinomial distribution
rho_ij = @(p_i,p_j) - p_i*p_j / sqrt(p_i*(1-p_i)*p_j*(1-p_j));

RhoA = zeros(length(p_a));
RhoB = zeros(length(p_b));

for i = 1:size(RhoA,1)
   for j = 1:size(RhoB,1)
       if (i ~= j)
           RhoA(i,j) = rho_ij(p_a(i),p_a(j));
           RhoB(i,j) = rho_ij(p_b(i),p_b(j));
       end
   end
end

figure;
subplot(1,2,1);
imagesc(RhoA); axis square; colorbar;

subplot(1,2,2);
imagesc(RhoB); axis square; colorbar;

figure;
imagesc(abs(RhoA-RhoB).^2); axis square; colorbar;

