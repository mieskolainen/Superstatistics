% Mobius/combinatorics inversion simulation
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src

% R to mu function
mu = @(R) -log(1 - R);

% Test
test = 'KS';   % Kolmogorov-Smirnov (default)
%test  = 'RMS'; % RMS
%test = 'KL';   % Kullback-Leibler

% Dirichtlet distribution concentration parameter 
% (-1 is for constant p = 1/n, which should actually give integrated
% results the same as 1.0 for Dirichlet)
dpam_vals = [0.5 1.0]; %[-1 0.1 1.0];

% Number of repetations of random simulations
Nsim = 100; % Order of 50-500

% Number of observables (dimensionality of the vector space)
N_vals   = [2 4 8]; % approx 14 is the maximum on standard PC

% Number of events per plot, units of 10^x
N_E_vals = [3 4 5];

% Mu-values range
mu_min = 1e-1;
mu_max = 30;
mu_steps = 12;
mu_values = logspace(log10(mu_min),log10(mu_max), mu_steps);

% Histogram values
MU_limit = 3;

% Fmincon options
options = optimoptions(@fmincon,'Display','iter','TolFun',1e-10,'TolX',1e-10);

% These are needed with optimization code
%global Nevents;
%global LAMBDA;
%global invLAMBDA;
%global C;
%global y;

% Different physical simulation scenarios
tic;
for dpam = dpam_vals

% Number of observables N = 2,5,8...
for w = 1:length(N_vals)
    
    % Open up figures
    f1 = figure;
    f2 = figure;
    
    % Legend texts of figure 1
    legends_f1 = {};
    legends_f2 = {};
    
    d = N_vals(w);       % Number of observables
    n = 2^d - 1;         % -> Number of components in binary vector

    % Create combinatorics matrix and remove ID 0
    C = createCBM(d);
    C = C(2:end,:);

    LAMBDA = sparse(amat(d));
    invLAMBDA = inv(LAMBDA);
    
    % Loop over trajectories of # events = 10^2, 10^4, 10^6 etc.
    for o = 1:length(N_E_vals)

        Nevents = 1*10^N_E_vals(o);
        
        % Simulation output
        X = zeros(2^d-1,Nsim,length(mu_values));
        Y = zeros(2^d-1,Nsim,length(mu_values));
        A = zeros(2^d-1,Nsim,length(mu_values));
        O = zeros(2^d-1,Nsim,length(mu_values));
    
        % Loop over mu-values
        for z = 1:length(mu_values)

            %global nu;
            nu = mu_values(z);           % Poisson mu-value
            
            % Loop over toy MC samples (PARALLEL LOOP here)
            parfor r = 1:Nsim
                
                fprintf('Simulation %d/%d, N = %d (dimension), E = %d [log10(#events)], alpha = %0.1f \n', r, Nsim, d, log10(Nevents), dpam);
                
                %% Create "physics model"
                
                % All components with equal probabilities
                if (dpam == -1)
                    x = ones(2^d-1,1); x = x / sum(x);
                    
                % Dirichlet distribution
                else
                    alpha = ones(1,2^d-1) * dpam;
                    x = dirnd(alpha,1)';  
                end

                % Make it column vector
                x = x(:);

                % Save input distribution
                X(:,z,r) = x;
                
                %% Create simulated data according to the fiducial rates
                MCtest_X = zeros(d,Nevents); % Allocate memory here
                MCtest_targets = randsample(1:2^d-1, Nevents, true, x)';
                for i = 1:Nevents
                    MCtest_X(:,i) = C(MCtest_targets(i),:)';
                end

                % Pileup simulation

                % Create the pileup mix
                tic
                [MC_pileup_X, MC_pileup_matrix, Nempty_BC] = poissonpileup_B(MCtest_X, MCtest_targets, nu, Nevents, size(C,1));

                % Make it binary
                MC_pileup_X = MC_pileup_X & ones(size(MC_pileup_X));

                % Classify the events
                MC_targets_pu = mybi2de(MC_pileup_X, 'right-msb');

                fprintf('Poisson \\mu-value: %0.3f \n', nu);
                %R_OR = (Nevents/(Nempty_BC + Nevents));
                %R_OR = (1-exp(-nu))*(1 + sqrt(Nevents)*randn(1)/Nevents); % Add noise
                R_OR = 1-exp(-nu);
                fprintf('MC (trigger) rate: %0.5f (oracle truth: %0.5f)\n', R_OR, 1-exp(-nu));
                toc
                
                %% Create the corresponding "smearing" matrix for debug, visualization etc.
%{
                if (o == length(length(N_E_vals)) && z == length(mu_values) && r == 1)
                    W = createsmearing(MC_targets_pu, MCtest_targets, C);
                    
                    f3 = figure;
                    imagesc(log10(W));% title('Smearing matrix W');
                    axis square;
                    set(gca,'fontsize',12);
                    pdfcrop(0.7, 0.96);
                    cmd = sprintf('print -dpdf ../figs/Wmat/mixing_N%d_mu%0.2f_alpha%0.1f.pdf', d, nu, dpam);
                    eval(cmd);
                    close(f3);
                end
%}
                
                %% Count the measurement data
                y = zeros(size(x));
                for i = 1:length(y)
                    y(i) = sum(MC_targets_pu == i);
                end
                y = y / sum(y);
                
                %% Error before correction

                e_before = norm(y - x);
                fprintf('Before inversion: L2-error: %0.6f \n', e_before);

                % Save output
                Y(:,z,r) = y;
                
                %% Analytic "all-order" pileup correction
                
                x_hat = analytic(y, R_OR, LAMBDA, invLAMBDA);
                
                % For a large Poisson
                % \mu we might get negative (non-sane) results => abs and renorm.
                x_hat = abs(x_hat); x_hat = x_hat / sum(x_hat);

                %% Error after inverse
                e_anal = norm(x_hat - x);
                fprintf('Analytic inversion: L2-error: %0.6f \n', e_anal);
                
                % Save result
                A(:,z,r) = x_hat;
                
                %% Solve as a fixed order optimization problem
                %{

                global loops;

                order = 5;

                % Create pileup combination matrices
                loops = cell(order,1);
                for i = 1:order
                    loops{i} = loop(2^d-1, i); 
                end

                x0 = x_hat;

                Aeq = ones(1,length(x0)); beq = 1; % Unitarity constraint, sum of prob = 1
                lb = zeros(size(x0)); % lower bound
                ub = ones(size(x0)); % upper bound

                x_hat_opt = fmincon(@nonlinear, x0, [], [], Aeq, beq, lb, ub, [],options);

                e_optim = norm(x_hat_opt - x);
                fprintf('Non-linear inversion: L2-error: %0.6f \n', e_optim);

                % Save results
                O(:,r,z) = x_hat_opt;

                x_hat_opt - x_hat

                fprintf('\n');
                %}

            end %% MC simulation

        end %% Mu-values

        % Save results
        %save(sprintf('pileupsim_%0.0f_d%d_NBC%d.mat', now, d, log10(Nevents)), 'X','Y','A','O','d','mu_values','Nevents');

        % Plotting data saved here
        m_A  = zeros(length(mu_values), 1);
        cl_A = zeros(length(mu_values), 2);
        
        m_Y  = zeros(length(mu_values), 1);
        cl_Y = zeros(length(mu_values), 2);

        % Histogram variables
        e_A_D = []; % Raw differences per component after correction
        e_Y_D = []; % Raw differences per component before correction
        
        % Plot figures
        for z = 1:length(mu_values)

            e_Y = zeros(Nsim,1); % Errors for before correction
            e_A = zeros(Nsim,1); % Errors for after correction
            
            X_ = squeeze(X(:,z,:));
            Y_ = squeeze(Y(:,z,:));
            A_ = squeeze(A(:,z,:));
            
            % Take histogram only for low mu-values, because at high mu we
            % are saturated -> no inversion
            if (mu_values(z) < MU_limit)
                % Calculate raw differences
                % ESTIMATE - TRUE or \hat{p} - p
                for r = 1:Nsim
                    
                    residual_Y = Y_(:,r) - X_(:,r);
                    residual_A = A_(:,r) - X_(:,r);
                    
                    % Save residuals
                    e_Y_D  =  [e_Y_D; residual_Y];
                    e_A_D  =  [e_A_D; residual_A];

                    % DEBUG
                    %fprintf('Before: %0.5f, After = %0.5f \n', norm(residual_Y), norm(residual_A));
                end
            end
            
            % Metrics
            if (strcmp(test,'KS')) % Kolmogorov-Smirnov
                for r = 1:Nsim
                    e_A(r) = max(abs(cumsum(X_(:,r)) - cumsum(A_(:,r))));
                    e_Y(r) = max(abs(cumsum(X_(:,r)) - cumsum(Y_(:,r))));
                end
            elseif (strcmp(test,'RMS')) % root mean square
                for r = 1:Nsim
                    e_A(r) = sqrt(mean((X_(:,r) - A_(:,r)).^2));
                    e_Y(r) = sqrt(mean((X_(:,r) - Y_(:,r)).^2));
                end
            elseif (strcmp(test,'KL')) % Kullback-Leibler
               for r = 1:Nsim
                    kl = A_(:,r).*log2(A_(:,r)./X_(:,r)); kl(isnan(kl) | isinf(kl)) = 0;
                    e_A(r) = sum(kl);
                    
                    kl = Y_(:,r).*log2(Y_(:,r)./X_(:,r)); kl(isnan(kl) | isinf(kl)) = 0;
                    e_Y(r) = sum(kl);
               end
            end
            
            % Median and percentiles of the "metric"
            m_A(z)    = median(e_A);
            cl_A(z,1) = m_A(z) - prctile(e_A,16);
            cl_A(z,2) = prctile(e_A,84) - m_A(z);

            m_Y(z)    = median(e_Y);
            cl_Y(z,1) = m_Y(z) - prctile(e_Y,16);
            cl_Y(z,2) = prctile(e_Y,84) - m_Y(z);
        end

        % To the same plot
        colors = {'k','b','r','g','y'};

        % HISTOGRAMS HERE
        BinEdges = linspace(-0.1,0.1,50);
        BinCenter = (BinEdges(1:end-1) + BinEdges(2:end))/2;
        fx = figure;
        hY = histogram(e_Y_D, BinEdges);
        hYvals = hY.Values;
        close(fx);
        fx = figure;
        hA = histogram(e_A_D, BinEdges);
        hAvals = hA.Values;
        close(fx);
        
        figure(f1);
        yvals = hYvals / sum(hYvals); % Normalize to discrete probability density
        avals = hAvals / sum(hAvals);
        stephist(BinCenter, yvals, [colors{o} ':']); hold on; % measurement Y
        stephist(BinCenter, avals, [colors{o} '-']); hold on;  % inverted A
        set(gca,'yscale','log');
        
        % MSE (mean squared error) in the legend
        legends_f1{end+1} = sprintf('$%0.1E$', mean(e_Y_D.^2));
        legends_f1{end+1} = sprintf('$%0.1E$', mean(e_A_D.^2)); 
        
        xlabel('$r_c$', 'interpreter', 'latex');
        ylabel('$f(r_c)$', 'interpreter','latex'); % TURNED OFF TO SAVE SPACE
        
        if (dpam ~= -1)
            title(sprintf('$N = %d, \\alpha = %0.1f, \\mu \\in [%0.1f,%0.0f]$', d, dpam, mu_values(1), MU_limit),'interpreter','latex');
        else
            % constant case
            title(sprintf('$N = %d, \\mu \\in [%0.1f,%0.0f]$', d, mu_values(1), MU_limit),'interpreter','latex');
        end
        
        axis square;
        xlim([-0.1 0.1]);              % manual scale x-axis
        ylim([prctile(avals, 1) 1]);   % manual scale y-axis
        
        % ----------------------------------------------------------------
        % CURVES HERE
        figure(f2);
        
        % Curves
        %x_vals = linspace(mu_values(1), mu_values(end), 5);
        %y_vals = ones(5,1)* (1/sqrt(Nevents)) * (1/sqrt(n));

        % Theoretical lower bound
        %loglog(x_vals, y_vals, 'color',[0.5 0.5 0.5],'linewidth', 1.20); hold on;

        % Values
        loglog(mu_values, m_Y(:), [colors{o} '.:']); hold on;
        loglog(mu_values, m_A(:), [colors{o} 's-']);
        
        %l = legend('Mixed','Inverted');
        %set(l,'location','northwest');
        %set(l,'Interpreter','Latex');
        %legend('boxoff');
        
        % Errorbars
        errorbar(mu_values, m_Y(:), cl_Y(:,1), cl_Y(:,2), [colors{o} '.:']);
        errorbar(mu_values, m_A(:), cl_A(:,1), cl_A(:,2), [colors{o} 's-']);

        xlabel('$\mu$','interpreter','latex');
        ylabel(sprintf('$%s$', test),'interpreter','latex');
        
        if (dpam ~= -1)
            title(sprintf('$N = %d, \\alpha = %0.1f$', d, dpam),'interpreter','latex');
        else % Const case
            title(sprintf('$N = %d$', d),'interpreter','latex');
        end
        axis square;
        if (~strcmp(test, 'KL'))
            axis([min(mu_values) max(mu_values) 1e-3 1]);
        else
            axis([min(mu_values) max(mu_values) 1e-5 100]);
        end
        
        %set(gca, 'XlimMode', 'auto');
        %set(gca,'XTick', mu_values)
        %grid on;        
    end %% N_E
    
    % Plot histograms
    figure(f1);
    l = legend(legends_f1);
    set(l,'Interpreter', 'Latex', 'location', 'northeast','fontsize',10);
    set(gca,'fontsize', 13);
    filename = sprintf('../figs/hist_N%d_alpha%0.1f.pdf', d, dpam);
    print(f1, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
    close(f1); % important
    
    % Plot curves
    figure(f2);
    set(gca,'fontsize', 13);
    filename = sprintf('../figs/%s_N%d_alpha%0.1f.pdf', test, d, dpam);
    print(f2, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));
    close(f2); % important
    
end % d-values
end % Dirichlet distribution scenarios

toc;
