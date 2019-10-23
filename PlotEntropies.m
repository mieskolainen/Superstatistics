% Plot entropy simulations
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;

addpath src

STEPS = 250; % Not too high, otherwise pdf turned to .png
mu_values = logspace(-1.5,log10(200),STEPS);

alpha_values = [0.01 0.1 0.2 1.0];

REALIZATIONS = 20; % Number of MC realizations
N_values = [2,4,8];

% mean and variance
f1 = {};
f2 = {};
for i = 1:length(alpha_values)
   f1{end+1} = figure;
   f2{end+1} = figure;
end

tic;
% Loop over dimensionality
for N = N_values
    
% Loop over dirichtlet distribution concentration parameters
for k = 1:length(alpha_values)

% Save entropy trajectories here
entropies = zeros(REALIZATIONS, STEPS);

% Loop over realizations
for a = 1:REALIZATIONS

% Dirichlet distribution input
alpha = ones(1,2^N-1) * alpha_values(k);
input = dirnd(alpha,1)';

% Create Möbius matrices
LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

y_values = zeros(length(input),STEPS);
p_values = zeros(length(input),STEPS);
entropy = zeros(1,STEPS);

% Loop over mu-values
for i = 1:length(mu_values)
    mu = mu_values(i);
    
    % Forward
    %
    y = LAMBDAINV*((exp(-mu*LAMBDA*input) - 1)/(exp(-mu)-1));
    y_values(:,i) = y;
    %}
    
    % Inverse
    %{
    p = LAMBDAINV*(log((exp(-mu) - 1)*LAMBDA*input + 1)) / (-mu);
    p_values(:,i) = p;
    %}
    entropies(a, i) = shannonentropy(y); % Shannon entropy in bits
end

end % Realizations loop

%
f3 = figure;
figure(f3);
plot(mu_values, entropies');

xlabel('$\mu$','interpreter','latex');
ylabel('$S(\mathbf{y})$','interpreter','latex');
title(sprintf('$N=%d, \\alpha = %0.2f$', N, alpha_values(k)), 'interpreter','latex');
set(gca,'xscale','log');
axis([min(mu_values) max(mu_values) 0 N]);
axis square;

    filename = sprintf('../figs/DIRrealizationsN%dk%d.pdf', N, k);
    print(f3, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
    
close(f3);
%}

% Calculate point statistics
figure(f1{k});
plot(mu_values, mean(entropies)); hold on;
set(gca,'xscale','log');
xlabel('$\mu$','interpreter','latex');
ylabel('$\langle S(\mathbf{y})\rangle$','interpreter','latex');
title(sprintf('$\\alpha = %0.2f$', alpha_values(k)), 'interpreter','latex');
axis square;
axis([min(mu_values) max(mu_values) 0 max(N_values)]);

figure(f2{k});
plot(mu_values, var(entropies)); hold on;
set(gca,'xscale','log');
xlabel('$\mu$','interpreter','latex');
ylabel('$\langle S(\mathbf{y})^2 \rangle - \langle S(\mathbf{y}) \rangle^2$','interpreter','latex');
title(sprintf('$\\alpha = %0.2f$', alpha_values(k)), 'interpreter','latex');
axis square;
axis([min(mu_values) max(mu_values) 0 inf]);

end % Alpha-values loop

end % N-loop
toc;

%% Make legends

legs = {};
for i = 1:length(N_values)
    legs{i} = sprintf('$N = %d$', N_values(i));
end

for k = 1:length(alpha_values)

figure(f1{k});
legend(legs,'interpreter','latex'); legend('boxoff');

    filename = sprintf('../figs/DIRmean%d.pdf', k);
    print(f1{k}, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

figure(f2{k});
legend(legs,'interpreter','latex'); legend('boxoff');
    
    filename = sprintf('../figs/DIRvar%d.pdf', k);
    print(f1{k}, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
    
end

