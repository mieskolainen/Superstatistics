% Cost function for rapidity gap data
%
% mikael.mieskolainen@cern.ch, 2019

function [totalchi2] = costfunc(x)

global param;
global setup;

% Setup parameters
param.npom       = x(1);

ptcuts    = [0.2, 0.4, 0.6, 0.8];

% Read in ATLAS data
DATA = cell(4,1);
cuts = cell(4,1);
parfor k = 1:4
    DATA{k}.table = readtable(sprintf('./HEPDATA/Table_%d.csv', k));
    
    % Set ATLAS cuts
    cuts{k}.etacut = 4.9;
    cuts{k}.ptcut = ptcuts(k);
end

chi2val = zeros(4,1);
obs = cell(4,1);

% Simulate over different pt cuts
parfor k = 1:4
    
    % Operate
    obs{k} = makesimulation(cuts{k}, param, setup);
    
    [~, chi2] = makehistogram(obs{k}, DATA{k}, param);
    %printhistogram();
    
    chi2val(k) = chi2;
end

totalchi2 = sum(chi2val);

fprintf('Parameters: \n');
disp(x);
fprintf('Total chi2 = %0.3f \n\n', totalchi2);

end