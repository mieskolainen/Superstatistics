% Gaps spacing statistics
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

mingap = @(n,alpha) (1-(n+1)*alpha).^n;

n = 1:10;
alphaval = linspace(1e-3,0.1,7);
legs = {};
for i = 1:length(alphaval)    
    f(:,i) = mingap(n, alphaval(i));
    legs{end+1} = sprintf('$\\alpha = %0.1E$', alphaval(i));
end
plot(n, f); hold on;
xticks(n);
l = legend(legs); set(l,'interpreter','latex');
axis square; axis tight; axis([min(n) max(n) 0 1]);
xlabel('$n$','interpreter','latex');
ylabel('$P(Y_j > \alpha)$','interpreter','latex');

%%

n = 1:10;
betaval = linspace(0.3,0.9,7);
legs = {};
for i = 1:length(betaval)
    f = maxgap(n, betaval(i));
    
    plot(n, f, 's-'); hold on;
    legs{end+1} = sprintf('$\\beta = %0.1f$', betaval(i));
end
xticks(n);
l = legend(legs); set(l,'interpreter','latex','location','southeast');
axis square; axis tight; axis([min(n) max(n) 0 1]);
xlabel('$n$','interpreter','latex');
ylabel('$P(Y_j < \beta)$','interpreter','latex');

