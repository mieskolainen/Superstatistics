% Distributions for different sampling statistics
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src

% Multinomial distribution (sampling with replacement)
m = @(p,x) factorial(sum(x))/prod(factorial(x)) * prod(p(:)'.^x);
%m = @(p,x) 1/prod(factorial(x)) * prod(p(:)'.^x);

% Entropy in bits
S = @(p) -sum(p .* log2(p));

N = 3; % Dimensionality

% 1. Random
%p = rand(1,7); p = p / sum(p);

% 2. Flat 
p = ones(1,7) + rand(1,7)*0.25; p = p / sum(p);

% 3. COMB input
%p = linspace(1/(2^N-1),1-1/(2^N-1),2^N-1)'; p = p'/sum(p); p = flipud(p);

% Poisson mean ticks
bins = 3000;
mu = linspace(0, 200, bins);

% Matrix for the "Bosonic" case
LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

for mode = 1:3

    Y_vals = zeros(7,size(mu,2));
    S_vals = zeros(1,size(mu,2));
    y = zeros(2^N-1,1);
    
    for i = 1:length(mu)

        if (mode == 1) % "Semi-Maxwell-Boltzmann"
            
            y = LAMBDAINV * ((exp(-mu(i)*LAMBDA*(p')) - 1)/(exp(-mu(i))-1));
            
        end
        if (mode == 2) % "Semi-Bosonic"
            
            % Pick probabilities (only three needed) from normalized Poisson distribution
            P = poisspdf([1:3], mu(i)); P = P / sum(P);
            
            y(1) = P(1)*p(1);
            y(2) = P(1)*p(2);
            y(3) = P(1)*p(3) + P(2)*(m(p,[1,1,0,0,0,0,0]) + m(p,[1,0,1,0,0,0,0]) + m(p,[0,1,1,0,0,0,0])) + P(3)*m(p,[1,1,1,0,0,0,0]);
            y(4) = P(1)*p(4);
            y(5) = P(1)*p(5) + P(2)*(m(p,[1,0,0,1,0,0,0]) + m(p,[1,0,0,0,1,0,0]) + m(p,[0,0,0,1,1,0,0])) + P(3)*m(p,[1,0,0,1,1,0,0]);
            y(6) = P(1)*p(6) + P(2)*(m(p,[0,1,0,1,0,0,0]) + m(p,[0,0,0,1,1,0,0]) + m(p,[0,1,0,0,1,0,0])) + P(3)*m(p,[0,1,0,1,0,1,0]);
            y(7) = 1 - sum(y(1:6));
        end
        if (mode == 3) % "Semi-Fermionic"
            
            % Pick probabilities (only three needed) from normalized Poisson distribution
            P = poisspdf([1:3], mu(i)); P = P / sum(P);
            
            y(1) = P(1)*p(1);
            y(2) = P(1)*p(2);
            y(3) = P(1)*p(3) + P(2)*m(p,[1,1,0,0,0,0,0]);
            y(4) = P(1)*p(4);
            y(5) = P(1)*p(5) + P(2)*m(p,[1,0,0,1,0,0,0]);
            y(6) = P(1)*p(6) + P(2)*m(p,[0,1,0,1,0,0,0]);
            y(7) = 1 - sum(y(1:6));
        end
        
        Y_vals(:,i) = y;
        S_vals(i)   = S(y); % Entropy
    end
    
    legends = {};
    for c = 1:length(y)
         legends{end+1} = sprintf('$y_{%d}$',c);
    end
    
    f1 = figure;
    plot(mu, Y_vals); hold on;
    %plot(mu, 1 - Y_vals(end,:), 'k-.'); % 1 - last
    
    xlabel('$\mu$','interpreter','latex');
    ylabel('$\mathbf{y}$','interpreter','latex');
    l = legend(legends); set(l,'interpreter','latex','location','southwest'); legend('boxoff');

    if (mode == 1)
        str = sprintf('Semi-MB');
    elseif (mode == 2)
        str = sprintf('Semi-BE');
    elseif (mode == 3) 
        str = sprintf('Semi-FD');
    elseif (mode == 4) 
        str = sprintf('AGK');
    end
    
    title(str,'interpreter','latex');
    
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    axis tight; axis square;
    axis([min(mu) max(mu) 1e-5 1]);
    xticks([0.1, 1, 10, 100]);
            
    %set(gca,'fontsize',13);
    
    filename = sprintf('../figs/statistics%d.pdf', mode);
    print(f1, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
    
    %close all;
end
