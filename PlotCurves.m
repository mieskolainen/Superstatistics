% Plot curves as a function of mu
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;

addpath src

N = 3;

% A.) Maximum Entropy input
%input = ones(2^N-1,1) + rand(2^N-1,1)*0.05; input = input/sum(input);

% B.) Comb [0...1]
input = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; input = input/sum(input);
input = flipud(input);

% B.) Comb [-1...1]
%input = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; 
%input = input - 0.5;
%input = flipud(input);

LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

STEPS = 1e3;
mu_values = logspace(-1.5,log10(30),STEPS);
y_values = zeros(length(input),STEPS);
p_values = zeros(length(input),STEPS);
entropy = zeros(1,STEPS);

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    % Forward
    y = LAMBDAINV*((exp(-mu*LAMBDA*input) - 1)/(exp(-mu)-1));
    y_values(:,i) = y;
    
    % Inverse
    p = LAMBDAINV*(log((exp(-mu) - 1)*LAMBDA*input + 1)) / (-mu);
    p_values(:,i) = p;
    
   % entropy(i) = -sum(y(:).*log(y(:)));
end

% loop over the to choose where y(end) ~ 0.65
target = 0.65;
min_ind = -1;
min_val = inf;
for j = 1:size(y_values,2)
    diff = abs( y_values(end,j) - target);
    if (diff < min_val)
        min_val = diff;
        min_ind = j;
    end
end
optim = y_values(:,min_ind);

f1 = figure;
plot(mu_values, y_values); hold on;
%plot(mu_values, entropy, 'k--');

xlabel('$\mu$','interpreter','latex');
ylabel('$y_c$','interpreter','latex');
%set(gca,'yscale','log');
set(gca,'xscale','log');
axis square; axis tight; 

legends = {};
for c = 1:length(y)
    legends{end+1} = sprintf('$y_{%d}$',c);
end
%legends{end+1} = sprintf('$H(y)$');

if (N <= 4)
    l = legend(legends);
    set(l,'interpreter','latex', 'location','northwest');
    legend('boxoff');
end
filename = sprintf('../figs/FLATp2y.pdf');
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

close all;

f2 = figure;
plot(mu_values, p_values); hold on;
plot(mu_values, zeros(length(mu_values),1), 'k-.'); % Horizontal axis
xlabel('$\mu$','interpreter','latex');
ylabel('$p_c$','interpreter','latex');
set(gca,'xscale','log');
%set(gca,'yscale','log');
axis tight; axis square;

legends = {};
for c = 1:length(y)
    legends{end+1} = sprintf('$p_{%d}$',c);
end

% Create smaller axes in top right, and plot on it
%{
% Restart color indexing
axes('Position',[0.28 0.60 .3 .3]); % x y
range = 0.9;
plot(mu_values(1:round(length(mu_values)*range)), p_values(:,(1:round(length(mu_values)*range))) );
axis([0,7, -0.015,0.04]); set (gca,'xscale','log');
set(gca,'ytick',[])

ax = gca;
ax.ColorOrderIndex = 1;
%}
%axis tight; box on;

if (N <= 4)
    l = legend(legends);
    set(l,'interpreter','latex', 'location','southwest');
    legend('boxoff');
end
filename = sprintf('../figs/FLATy2p.pdf');
print(f2, filename, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

