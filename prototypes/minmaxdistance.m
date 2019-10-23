% Minimum and Maximum gap spacing distributions / stick lengths
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

S        = 1e5;
n_values = [2 3 4 5 6];

% Number of particles
legs = {};
for n = n_values

    legs{end+1} = sprintf('$n = %d$', n);
    
    mind = zeros(S,1);
    mead = zeros(S,1);
    maxd = zeros(S,1);
    
    for ind = 1:S
        
        % Uniform particles over rapidity
        y = rand(n,1);
        
        % Order
        y = sort(y);
        
        % Calculate particle-to-particle distances
        delta = zeros(n-1, 1);
        k = 1;
        for i = 1:n-1
           delta(k) = abs(y(i) - y(i+1));
           k = k + 1;
        end
        
        % Calculate pairwise distances
        %{
        delta = zeros((n^2-n)/2, 1);
        k = 1;
        for i = 1:n
            for j = i+1:n
                delta(k) = abs(y(i) - y(j));
                k = k + 1;
            end
        end
        %}
        
        % Point statistics
        mind(ind) =  min(delta);
        mead(ind) = mean(delta);
        maxd(ind) =  max(delta);
    end
    
    xedge = linspace(0,1,100);
    
    subplot(1, 3, 1);
    c = hist1m(mind, xedge);
    stephistedge(xedge, c/sum(c), '-');
    axis square;
    xlabel('$\Delta_{min}$','interpreter','latex');
    ylabel('$P$','interpreter','latex');
    hold on;
    l = legend(legs); set(l,'interpreter','latex');
    
    subplot(1, 3, 2);
    c = hist1m(mead, xedge);
    stephistedge(xedge,c/sum(c), '-');
    axis square;
    xlabel('$\Delta_{mean}$','interpreter',   'latex');
    ylabel('$P$','interpreter','latex');
    hold on;

    subplot(1, 3, 3);
    c = hist1m(maxd, xedge);
    stephistedge(xedge,c/sum(c), '-');
    axis square;
    xlabel('$\Delta_{max}$','interpreter','latex');
    ylabel('$P$','interpreter','latex');
    hold on;
end


%% Analytic stick lengths

figure;

d = linspace(0,1,1e3);
Y = 1; % Interval

for n = 3:7 % Number of intervals
    
    P = sticklength(d,n);
    subplot(1,2,1);
    
    j = 1;
    
    plot(d, P(j,:)); hold on; axis square;
    xlabel('$\Delta_{max}$','interpreter','latex');
    ylabel('$P$','interpreter','latex'); axis tight;
    
    subplot(1,2,2);
    plot(d(1:end-1), diff(P(j,:))/(d(2)-d(1))); hold on; axis square;
    xlabel('$\Delta_{max}$','interpreter','latex');
    ylabel('$P$','interpreter','latex'); axis tight;
end



