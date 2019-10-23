%% ANALYTIC Distribution with p_c = 1/n (max input entropy) constraint target
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;

addpath src

N = 3;

mu_values = logspace(log10(1e-3),log10(30), 1e3);

f1 = figure;

y_values = zeros(7,length(mu_values));
g_values = zeros(1,length(mu_values));
h_values = zeros(1,length(mu_values));

n = 7;

for i = 1:length(mu_values)
    mu = mu_values(i);
    e1  = exp(mu)-1; % e_^mu

    y(1) = (exp(mu/n)-1)/e1;
    y(2) = y(1);
    y(3) = (-2*exp(mu/n) + exp(3*mu/n) + 1)/(exp(mu)-1);
    y(4) = y(1);
    y(5) = y(3);
    y(6) = y(3);
    y(7) = 1 - 3*y(1) - 3*y(3);
    y_values(:,i) = y;
    
    g_values(i) = sum(y([3 5 6])) - sum(y([1 2 4]));
    h_values(i) = sum(y(1:6));
end

plot(mu_values, y_values); hold on;
set(gca,'xscale','log');
axis tight; axis square;
xlabel('$\mu$','interpreter','latex');
ylabel('$\mathbf{y}$','interpreter','latex');

legends = {};
for c = 1:length(y)
     legends{end+1} = sprintf('$y_{%d}$',c);
end
legends{end+1} = sprintf('$y_{3,5,6} - y_{1,2,4}$'); % g
legends{end+1} = sprintf('$1 - y_7$'); % h

% Sum distributions
range = 1;
plot(mu_values(1:range*length(mu_values)), g_values(1:range*length(mu_values)), 'k-'); hold on;
plot(mu_values(1:range*length(mu_values)), h_values(1:range*length(mu_values)), 'k:'); hold on;

l = legend(legends);
%set(l,...
%    'Position',[0.615178571428572 0.392222223584615 0.18422615414574 0.327619046256656],...
%    'Interpreter','latex');
set(l,'location','northwest','interpreter','latex');
legend('boxoff');

filename = sprintf('../figs/maxinent.pdf');
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

