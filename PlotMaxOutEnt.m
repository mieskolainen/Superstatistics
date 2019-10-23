% Maximum output entropy: Plot curves as a function of Poisson mu
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src

for N = 5:8

% A.) Maximum Entropy input
input = ones(2^N-1,1) + rand(2^N-1,1)*0.01; input = input/sum(input);

% B.) Comb [0...1]
%input = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; input = input/sum(input);
%input = flipud(input);

% B.) Comb [-1...1]
%input = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; 
%input = input - 0.5;
%input = flipud(input);

LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

STEPS = 200; % Not too high, otherwise pdf turned to scalar picture
mu_values = logspace(-1.5,log10(10),STEPS);
y_values = zeros(length(input),STEPS);
p_values = zeros(length(input),STEPS);
entropy = zeros(1,STEPS);

for i = 1:length(mu_values)
    mu = mu_values(i);
    
    % Forward
    %{
    y = LAMBDAINV*((exp(-mu*LAMBDA*input) - 1)/(exp(-mu)-1));
    y_values(:,i) = y;
    %}
    
    % Inverse
    p = LAMBDAINV*(log((exp(-mu) - 1)*LAMBDA*input + 1)) / (-mu);
    p_values(:,i) = p;
    
   % entropy(i) = -sum(y(:).*log(y(:)));
end

f1 = figure;

% Rectangle
h0 = rectangle('Position',[0 -1 max(mu_values) 1.0], 'FaceColor', ones(3,1)*0.96, 'Edgecolor', ones(3,1)); hold on;

h1 = plot(mu_values, zeros(length(mu_values),1), 'k-.'); hold on; % Horizontal axis
h2 = plot(mu_values, p_values);

xlabel('$\mu$','interpreter','latex');
ylabel('$\mathbf{p}$','interpreter','latex');

%set(gca,'yscale','log');
%set(gca,'xscale','log');
axis tight; axis square;
axis([-inf inf -0.05 0.1]);
xticks([0:1:max(mu_values)]);
yticks([linspace(-0.05, 0.1, 7)]);

set(gca, 'Layer', 'top'); grid off; % axis on top
set(gca,'box','on');

legends = {};
for c = 1:length(p)
    legends{end+1} = sprintf('$p_{%d}$',c);
end
title(sprintf('$N = %d$', N),'interpreter','latex');

if (N <= 3)
    l = legend(h2, legends);
    set(l,'interpreter','latex', 'location','southwest');
    legend('boxoff');
end

% Create smaller axes in top right, and plot on it
%{
if (N > 4)
    % Restart color indexing
    axes('Position',[0.241403508771931 0.595641646489104 0.342807017543859 0.306779661016956]); % x y
    range = 0.9;
    plot(mu_values(1:round(length(mu_values)*range)), p_values(:,(1:round(length(mu_values)*range))) ); hold on;
    plot(mu_values, zeros(length(mu_values),1), 'k-.'); % Horizontal axis
    axis([0,7, -0.015,0.04]); set (gca,'xscale','log');
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);

    ax = gca;
    ax.ColorOrderIndex = 1;
end
%}
%axis tight; box on;

filename = sprintf('../figs/maxoutentN%d.pdf', N);
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
close all;

end

