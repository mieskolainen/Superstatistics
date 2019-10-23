% Distributions for the negative ('unphysical') Poisson mu case
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src

% Multinomial distribution (sampling with replacement)
m = @(p,x) factorial(sum(x))/prod(factorial(x)) * prod(p(:)'.^x);

% Entropy in bits
S = @(p) -sum(p .* log2(p));

% Poisson mean ticks (not too many, otherwise pdf -> scalar png)
bins = 350;
mu = linspace(-100, 30, bins);

for N = 2:8 % Dimensionality

% 1. Random
%p = rand(1,2^N-1); p = p / sum(p);

% 2. Flat with tiny noise to distinquish between different modes
input = ones(1,2^N-1) + 0.05*randn(1,2^N-1); input = input / sum(input);

% 3. COMB input
%p = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; p = p'/sum(p); p = flipud(p);


% Matrix for the "Bosonic" case
LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

for mode = 1

    Y_vals = zeros(2^N-1,size(mu,2));
    S_vals = zeros(1,size(mu,2));
    y = zeros(2^N-1,1);
    
    for i = 1:length(mu)
        
        y = LAMBDAINV * ((exp(-mu(i)*LAMBDA*(input')) - 1)/(exp(-mu(i))-1));        
        %p = LAMBDAINV*(log((exp(-mu(i)) - 1)*LAMBDA*(input') + 1)) / (-mu(i));
        
        Y_vals(:,i) = y;
        S_vals(i)   = S(y); % Entropy
    end
    
    legends = {};
    for c = 1:length(y)
         legends{end+1} = sprintf('$y_{%d}$',c);
    end
    
    f1 = figure;
    plot(mu, Y_vals); hold on;
    plot(mu, 1 - sum(Y_vals,1), 'k-.');
    
    xlabel('$\mu$','interpreter','latex');
    ylabel('$\mathbf{y}$','interpreter','latex');
    
    if (N <= 3)
    l = legend(legends); set(l,'interpreter','latex','location','southwest'); legend('boxoff');    
    end
    
    %set(gca,'yscale','log');
    %set(gca,'xscale','log');
    axis tight; axis square;
    axis([min(mu) max(mu) -1.1 1.1]);
    set(gca,'XTick',min(mu):20:max(mu));
    title(sprintf('$N = %d$', N),'interpreter','latex');
    %set(gca,'fontsize',13);
    
    %set(gca,'yscale','log');
    
    filename = sprintf('../figs/negativemuN%d.pdf', N);
    print(f1, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
    
    close all;
end
end

