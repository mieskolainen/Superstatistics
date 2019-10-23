% Synthetic (toy) simulation of increasing pT and N_ch (multiplicity)
% threshold for partial cross sections with uniform dN/dy
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src
addpath /home/user/cernbox/#matlabcodes/

% Number of rapidity slices
% Construct algebraic alternating sign matrix
N = 3;
M_ = mmat(N);

% Parameters
meanpt = 1; % GeV
s0 = meanpt / (sqrt(pi/2));

% Number of MC events
N_MC = 1e4;
sigma_tot = 1;

% pT and N (multiplicity) thresholds
steps = 100;
pt_threshold = logspace(-1, 1, steps);
N_threshold  = logspace(-1, log10(40), steps);

% Poisson mean multiplicity per rapidity interval
lambdas = [0.5 0.8 1.0 2.0 4.0];

% DATA
PDATA = zeros(length(lambdas), 2^N, length(pt_threshold));
NDATA = zeros(length(lambdas), 2^N, length(N_threshold));

% Loop over different Poisson means
for w = 1:length(lambdas)
    
    % Threshold loop in PARALLEL
    parfor z = 1:steps
        
        fprintf('Lambda = %d/%d, Threshold %d/%d \n', w, length(lambdas), z, steps);
        sigma_pt = zeros(2^N,1);
        sigma_N  = zeros(2^N,1);
        
        % Event loop
        for i = 1:N_MC

            % Create toy transverse momentum (pt) measurement, by drawing
            % random number of charged particles ber rapidity interval (random
            % variable X)
            
            % Bernoulli random variables
            I_pt = zeros(N,1);
            I_N  = zeros(N,1);
            
            % Loop over rapidity slices, make measurements
            for k = 1:N
            % Create synthetic particles from Rayleigh pT distribution
                X = raylrnd(s0, 1, poissrnd(lambdas(w), 1));

            % Apply detector pT threshold
                T_pt = heaviside(X - pt_threshold(z));

            % Apply detector N (multiplicity) threshold
                T_N = heaviside(length(X) - N_threshold(z));

            % Sum over and turn into Bernoulli (binary) random variable
                I_pt(k) = sum(T_pt,2) & 1;
                I_N(k)  = sum(T_N,2) & 1;
            end
    
            % ------------------------------------------------------------
            % Algebraic construction
            
            % Kronecker products (1 X_1)^T (*) (1 X_2)^T (*) ... (*) (1 X_N)^T 
            K_pt = kron([1; I_pt(1)],[1; I_pt(2)]);
            K_N  = kron([1; I_N(1)], [1; I_N(2)]);
            for k = 3:N
               K_pt = kron(K_pt, [1; I_pt(k)]);
               K_N  = kron(K_N,  [1; I_N(k)]);
            end
            
            % Map with M_ matrix
            sigma_pt = sigma_pt + sigma_tot * M_* K_pt;
            sigma_N  = sigma_N  + sigma_tot * M_* K_N;
            % ------------------------------------------------------------    
        end
        
        PDATA(w,:,z) = sigma_pt / N_MC; % MC mean values
        NDATA(w,:,z) = sigma_N  / N_MC;
    end
end


%% PLOT

% Create legends
[legends, mark] = makelegend(N);

FONT = 13;

for w = 1:length(lambdas)
    for m = 1:2
        figure;
        for i = 1:2^N
            if (m == 1)
                stephist(pt_threshold / meanpt, PDATA(w,i,:), mark{i},'linewidth',1.1); hold on;
            elseif (m == 2)
               stephist(N_threshold / lambdas(w), NDATA(w,i,:), mark{i},'linewidth',1.1); hold on; 
            end
        end
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        ylabel('probability','interpreter','latex','fontsize',FONT); axis tight;
        title(sprintf('$\\langle N_{ch}/\\Delta_{\\eta} \\rangle = %0.1f$', lambdas(w)),'interpreter','latex');
        l = legend(legends);
        set(l,'interpreter','latex','fontsize',9); legend('boxoff');
        
        filename = '';
        if (m == 1)
            xlabel('$z = p_t^{cut} / \langle p_t \rangle$','interpreter','latex','fontsize',FONT);
            filename = sprintf('../figs/toysim_pt_%d.pdf', w);
        elseif (m == 2)
            xlabel('$z = N_{ch}^{cut} / \langle N_{ch} \rangle$','interpreter','latex','fontsize',FONT);
            filename = sprintf('../figs/toysim_N_%d.pdf', w);
        end
        %axis square;
        
        cmd = sprintf('print -dpdf %s', filename);
        eval(cmd);
        system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
        
        pdfcrop(0.93,1.01);
        eval(cmd);
        close all;
    end
end

