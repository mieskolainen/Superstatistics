% Test relative Entropy scaling for two dirichlet distributions with
% respect to the flat distribution
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

d_vals = [1:12];
N_MC = 5000;

% KL-divergence (relative entropy) base-2
KL_pq = @(p,q) sum(p.*log2(p./q));

% Alpha-values
alpha_vals = [0.1 1];

% Values
klvals  = zeros(length(alpha_vals), length(d_vals), N_MC);
klvals2  = zeros(length(alpha_vals), length(d_vals), N_MC);


for i = 1:length(d_vals) % Over dimensions
    
    d = d_vals(i);
    
    fprintf('Simulation %d/%d (d = %d) \n', i, length(d_vals), d);
    
    % Flat distribution
    p = ones(2^d-1,1)/(2^d-1);
    
    % Loop MC runs
    for k = 1:N_MC
        
        for a = 1:length(alpha_vals)
            
            % Random distribution from Dirichlet
            alpha = ones(1,2^d-1) * alpha_vals(a);
            q = dirnd(alpha,1)'; 

            kl = KL_pq(p, q); kl(isnan(kl) | isinf(kl)) = 0;
            klvals(a,i,k) = kl;
           
            % Reverse p and q
            kl = KL_pq(q, p); kl(isnan(kl) | isinf(kl)) = 0;
            klvals2(a,i,k) = kl;
        end
    end
end

%%
figure;
for i = 1:length(alpha_vals)
    
    subplot(1,2,i);
    
    for z = 1:1

        if (z == 1)
            values = squeeze(klvals(i,:,:));
            marker = '-';
        else
            values = squeeze(klvals2(i,:,:));
            marker = '--';
        end
        
        errorbar(d_vals, median(values,2), ...
            median(values,2)' - prctile(values', 95), ...
            median(values,2)' - prctile(values', 5), ['r' marker]);
        hold on;
        errorbar(d_vals, median(values,2), ...
            median(values,2)' - prctile(values', 84), ...
            median(values,2)' - prctile(values', 16), ['ks' marker]);

        axis square; axis tight;
        hold on;
        
        xlabel('$N$','interpreter','latex');
        if ( z==1)
            ylabel('KL$(p|q)$','interpreter','latex');
        else
            ylabel('KL$(q|p)$','interpreter','latex');
        end
        title(sprintf('$p_c = 1/n, q \\sim D(\\alpha = %0.1f)$', ...
        alpha_vals(i)),'interpreter','latex');
    end
    
    % Calculate entropies of the flat distribution
%     H = zeros(length(d_vals),1);
%     for j = 1:length(d_vals)
%         p = ones(2^d_vals(j)-1,1); p = p / sum(p);
%         H(j) = -sum(p.*log2(p));
%     end
    % plot(d_vals, H, 'r--');
end

%pdfcrop(0.9, 0.6);
%eval(sprintf('print -dpdf ../figs/KL_test.pdf'));
