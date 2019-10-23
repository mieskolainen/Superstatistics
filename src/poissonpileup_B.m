% Poisson distributed pile-up generator (version B)
%
% INPUT:        X   =   D data vectors (d dim) for K scattering processes,
%                       mixed with cross-sections, d x N
%        X_targets  =   Classification (N x 1)
%           lambda  =   Poisson distribution mean value, i.e. pile-up
%              N_Y  =   Number of vectors to be generated
%
% OUTPUT:        Y  =   Mixed data
%            Y_target_matrix
%
% mikael.mieskolainen@cern.ch, 2019

function [Y, Y_target_matrix, Nempty_BC] = poissonpileup_B(X, X_targets, nu, N_Y, K)

[d,N_X] = size(X);
Y = zeros(d,N_Y);
Y_target_matrix = zeros(N_Y,K);

% ------------------------------------------------------------------------
% Create Poisson probabilities of k = 1,2,...,k_max for fast random numbers
max_k = 300;
poiss_pdf = poisspdf(1:max_k, nu);
poiss_pdf = poiss_pdf / sum(poiss_pdf);
poissons = randsample(1:max_k, N_Y, true, poiss_pdf);
% ------------------------------------------------------------------------

Nempty_BC = 0; % True empty beam crossings (with exact Poisson randoms)

% Event loop
for i = 1:N_Y
    
% Exact Poisson random numbers, very slow
%     % Draw number of inelastic events for this bunch cross
%     n = poissrnd(nu);
%     
%     % Interaction condition
%     while (n < 1)
%         Nempty_BC = Nempty_BC + 1; % Add one empty-empty beam crossing
%         n = poissrnd(nu);
%     end
    
    % First put the original
    Y(:,i) = X(:,i);
    Y_target_matrix(i, X_targets(i)) = 1;
    
    % If we have pile-up, poissons(i) is the number of pileup
    if ( poissons(i) >= 2 )
        
        % Draw event classes, n-1 because we already chose one
        ind = randi(N_X, [poissons(i)-1 1]);
        
        % Vector
        Y(:,i) = Y(:,i) + sum(X(:,ind),2);
        
        % Targets
        for j = 1:poissons(i) - 1
            Y_target_matrix(i, X_targets(ind(j))) = ...
                Y_target_matrix(i, X_targets(ind(j))) + 1;
        end
    end
    
end
    
end
%}