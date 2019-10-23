% Good-Walker Toy Model style Mueller
%
% See G. Gustafson,
%    "The relation between the Good-Walker
%     and triple-regge formalism for diffractive excitation"
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

% Rapidity range
y = linspace(0,30,1e2);

DELTA = 0.2; % Pomeron intercept - 1
n = 1:10;    % Discrete multiplicities

% The solution to the toy evolution
P = zeros(length(n), length(y));
for i = 1:length(n)
   P(i,:) = exp(-DELTA*y).*(1-exp(-DELTA*y)).^(n(i)-1);
end

legends = {};
for c = 1:length(n)
     legends{end+1} = sprintf('$n = %d$',c);
end

plot(y, P);
set(gca,'xscale','log');
set(gca,'yscale','log');

%l = legend(legends); set(l, 'interpreter','latex'); legend('boxoff');

xlabel('$y$', 'interpreter','latex');
ylabel('$P_n(y)$', 'interpreter','latex');
axis square; axis tight;

%%
% Same on Bjorken-x, y = ln(1/x) <=> x = exp(-y)
% The parton distribution as the average number of partons at given x
% Bjorken
% xG(x) = <n> = \sum_n n P_n(y)

xGx = @(x, alpha) (1./x).^(DELTA);
x = logspace(-8,0,1e4);set(gca,'xscale','log');

figure;
plot(x, xGx(x, DELTA));
xlabel('$x$','interpreter','latex');
ylabel('$xG(x)$','interpreter','latex');
set(gca,'xscale','log');
%set(gca,'yscale','log');
axis square;

