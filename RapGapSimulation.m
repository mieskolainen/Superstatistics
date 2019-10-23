% Superfast Rapidity Gap Distribution Simulation with SD,DD,ND
%
% mikael.mieskolainen@cern.ch, 2019

clear; % Use this for globals!
close all;

addpath src
addpath /home/user/cernbox/#matlabcodes

%%
global setup;
global param;

% Simulation setup
setup.NMC          = 1e6;           % Number of accepted events for simulation
setup.mandelstam_s = 7000^2;        % CMS energy squared (GeV^2)

% Free (dynamical) paramaters
param.deltaP       = 0.09;          % Pomeron intercept is: 1 + deltaP
param.npom         = 5.0;           % Particle density (non-diffractive)
param.diff         = [0.27, 0.15];  % Diffraction 2xSD, DD
param.sigmainel    = 75;            % Total inelastic

param.lowmass      = 1.2;           % Minimum diffractive mass (GeV)
param.XI_max       = 1.0;           % Maximum diffractive mass

% FIT
%{
x0 = [param.npom];
xhat = fminsearch(@costfunc, x0);
%}

%% SIMULATE AND HISTOGRAM

ptcuts = [0.2, 0.4, 0.6, 0.8];

% Read in data
DATA = cell(4,1);
cuts = cell(4,1);
for k = 1:4
    DATA{k}.table = readtable(sprintf('./HEPDATA/Table_%d.csv', k));
    
    % Set cuts
    cuts{k}.etacut = 4.9;
    cuts{k}.ptcut = ptcuts(k);
end

chi2val = zeros(4,1);
obs     = cell(4,1);
histo   = cell(4,1);

parfor k = 1:4
    
    % Operate
    obs{k} = makesimulation(cuts{k}, param, setup);
    [histo{k}, chi2val(k)] = makehistogram(obs{k}, DATA{k}, param);
end
fprintf('Total chi2 = %0.3f \n', sum(chi2val));

%% Print
for k = 1:4
    printhistogram(histo{k}, param, cuts{k}) 
end
close all;
