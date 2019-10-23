% Create rapidity gap simulations
%
% mikael.mieskolainen@cern.ch, 2018

function obs = makesimulation(cuts, param, setup)

% SIMULATION

% Mean multiplicity ansatz given system with mass squared M2
diffractive_N    = @(M2) (1 + log(M2));

% Proton mass
mp     = 0.938;
etaMIN = -log(sqrt(setup.mandelstam_s)/mp);
etaMAX = -etaMIN;

% Observables
obs.observable    = cell(2,1);
obs.observable{1} = zeros(setup.NMC, 2);
obs.observable{2} = zeros(setup.NMC, 2);

% Get relative probabilities for process selection
frac = [1-sum(param.diff), param.diff] / param.sigmainel; % fractions of ND, 2XSD, DD
trials   = 0;
accepted = 0;

% Predraw random numbers for increased efficiency
poissonrandom = poissrnd(param.npom, setup.NMC*2, 1);
poissonrandom_k = 1;

% Monte Carlo Loop
while (true)
    
    % Pick process type
    while (true)
        ID = randi(3);
        if (rand(1) < frac(ID))
            break;
        end
    end
    
    %% Non-Diffractive
    if (ID == 1) 
      
    % Draw number of exchanges
    npom = poissonrandom(poissonrandom_k); poissonrandom_k = poissonrandom_k + 1;
    npom = max(1,npom);
    
    n = 0;
    for i = 1:npom
        NN = log(setup.mandelstam_s);
        nthis = round(poissrnd(NN, 1, 1));
        n = n + nthis;
    end
    
    %{
    % Draw number of exchanges
    npom = poissrnd(param.npom*(etaMAX - etaMIN), 1, 1);
    npom = max(1,npom);
    
    % Sample from lognormal
    n = sum(logrnd(param.loggamma, npom));
    %}
    
    % Sort particle over pseudorapidity
    etas = UNIF(etaMIN, etaMAX, n);
    etas = sort(etas);
    
    end
    
    %% SD (Single Diffractive) process simulation
    if (ID == 2)
        
    M2 = POWERF(param.lowmass^2, param.XI_max * setup.mandelstam_s, -(1+param.deltaP), 1); % Pick invariant mass squared
    rapgap = -log(M2 / setup.mandelstam_s);                     % Corresponding average rapidity gap
    
    % Pick diffractive system multiplicity
    N = diffractive_N(M2);
    n_X = max(1, poissrnd(N)); % Require at least one charged final state
    
    etas = UNIF(etaMIN, (etaMAX - rapgap) , n_X);
    etas = sort(etas);
    
    % Add one forward proton +rapidities
    etas = [etas; etaMAX];
    
    % Choose forward/backward single diffraction (parity mirror)
    if (rand(1) < 0.5)
        etas = -etas;
        etas = sort(etas);
    end
    
    end
    
    %% DD (Double Diffractive) process simulation
    if (ID == 3)
    
    while (true)
        
        r1 = POWERF(param.lowmass^2, setup.mandelstam_s, -(1+param.deltaP), 1);
        DD_M2_max = param.XI_max * (mp^2 * setup.mandelstam_s) / r1;
        r2 = POWERF(param.lowmass^2, DD_M2_max, -(1+param.deltaP), 1);
        
        % Random permutation
        if (rand(1) < 0.5)
            M2_X = r1;
            M2_Y = r2;
        else
            M2_Y = r1;
            M2_X = r2;
        end
    	% Check do we obey the limit in total
        if ((M2_X * M2_Y) / (mp^2 * setup.mandelstam_s) < param.XI_max)
            break;
        end
    end
        
    rapgap_X = -log(M2_X / setup.mandelstam_s);                   % Corresponding average rapidity gap
    rapgap_Y = -log(M2_Y / setup.mandelstam_s);                   % Corresponding average rapidity gap
    
    % Pick diffractive system multiplicity
    N = diffractive_N(M2_X);
    n_X = max(1, poissrnd(N)); % Require at least one charged final state
    
    % Pick diffractive system multiplicity
    N = diffractive_N(M2_Y);
    n_Y = max(1, poissrnd(N)); % Require at least one charged final state
    
    % Draw over rapidity
    etas_X = UNIF(etaMIN, (etaMAX - rapgap_X), n_X);   % left
    etas_Y = UNIF((etaMIN + rapgap_Y), etaMAX, n_Y);   % right
    etas   = [etas_X; etas_Y];
    etas   = sort(etas);
    
    end
    
    %% Generate transverse momentum
    
    pts = tsallisrnd(1.13, 0.084, 2.5, length(etas));
    
    %{
    % Brownian (Rayleigh) diffusion process
    if (rand(1) < 0.3 )
    gamma = 0.38;
    pts = raylrnd(gamma, [length(etas) 1]);
    else
    
    %
    % Cauchy superdiffusion process
    
    px  = CAUCHYRND(0, param.gamma, length(etas));
    py  = CAUCHYRND(0, param.gamma, length(etas));
    pts = sqrt(px.^2 + py.^2);
    end
    %}
    
    
    %% Apply fiducial cuts
    
    ind = (pts > cuts.ptcut) & (etas < cuts.etacut) & (etas > -cuts.etacut);
    
    etas = etas(ind);
    pts  = pts(ind);
    
    %% Event inside fiducial phase space
    
    if (~isempty(etas))
        
        % Calculate observables
        
        % 1. Forward
        delta_F = cuts.etacut - max(etas);
        obs.observable{1}(accepted+1,:) = [delta_F, ID];         
        
        % 2. Center (floating) gap
        all = [-cuts.etacut; etas; cuts.etacut];                    % Add boundaries
        %all = etas;
        delta = max( abs(all(2:end) - all(1:end-1)) ); 
        obs.observable{2}(accepted+1,:) = [delta, ID];
        
        % 3. Forward/Backward maximum
        %delta_B = abs(min(etas) - (-etacut));
        %delta   = max(delta_F, delta_B);
        %obs.observable{3}(accepted+1,:) = [delta, ID];
        
        accepted = accepted + 1;
    end
    
    %% Check event statistics
    
    trials = trials + 1;
    
    % Generated enough events
    if (accepted == setup.NMC)
        break;
    end
end

% Total fiducial inelastic
obs.fiducial_inel = param.sigmainel * (accepted / trials);
fprintf('ptcut = %0.1f GeV :: total fiducial inelastic = %0.1f mb \n', cuts.ptcut, obs.fiducial_inel);

obs.trials = trials;
obs.accepted = accepted;

end