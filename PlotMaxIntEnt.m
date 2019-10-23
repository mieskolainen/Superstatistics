% Maximum input entropy: Plot curves as a function of Poisson mu
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;

addpath src

STEPS = 200; % Not too high, otherwise pdf turned to scalar picture
mu_values = logspace(-1,log10(50),STEPS);
entropies = zeros(8, STEPS);

for N = 2:8

% A.) Maximum Entropy input
%input = ones(2^N-1,1) + rand(2^N-1,1)*0.0; input = input/sum(input);

% B.) Comb [0...1]
%input = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; input = input/sum(input);
%input = flipud(input);

% C.) Create basis vector input
C = createCBM(N);
C = C(2:end,:);
CSUM = sum(C,2);
ind = find(CSUM == 1);
input = zeros(2^N-1,1);
input(ind) = 1;
input = input / sum(input);

% Matrix
LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

y_values = zeros(length(input),STEPS);
p_values = zeros(length(input),STEPS);
entropy = zeros(1,STEPS);

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
    entropies(N, i) = shannonentropy(y);
end


%%

f1 = figure;
h1 = plot(mu_values, zeros(length(mu_values),1), 'k-.'); hold on; % Horizontal axis
h2 = plot(mu_values, y_values);

xlabel('$\mu$','interpreter','latex');
ylabel('$\mathbf{y}$','interpreter','latex');
set(gca,'xscale','log');
set(gca,'yscale','log');
axis tight; axis square;
axis([-inf inf 1e-3 1.0]);

legends = {};
for c = 1:length(y)
    legends{end+1} = sprintf('$y_{%d}$',c);
end
title(sprintf('$N = %d$', N),'interpreter','latex');

if (N <= 4)
    l = legend(h2, legends);
    set(l,'interpreter','latex', 'location','northeast');
    legend('boxoff');
end

% Create smaller axes in top right, and plot on it
%{
if (N > 4)
% Restart color indexing
axes('Position',[0.28 0.60 .3 .3]); % x y
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

filename = sprintf('../figs/maxinentN%d.pdf', N);
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
%close all;

end

%% Shannon entropy scaling (bits)

f2 = figure;

Nvals = 2:8;
plot(mu_values, entropies(Nvals,:)');
set(gca,'xscale','log');

legs = {};
for i = 1:length(Nvals)
   legs{i} = sprintf('$N=%d$',Nvals(i)); 
end

l = legend(legs); set(l,'interpreter','latex','location','northwest'); legend('boxoff');
xlabel('$\mu$','interpreter','latex');
ylabel('$S(\mathbf{y})$','interpreter','latex');
axis square;
axis([mu_values(1) mu_values(end) 0 9]);

filename = sprintf('../figs/entropyscaling.pdf');
print(f2, filename, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

