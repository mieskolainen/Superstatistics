% Number of edges in N-dimensional hypercube, oeis.org A001787
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

% Loop
counts = [];
N_values = 1:6;
for N = N_values
C = createCBM(N);
counts(end+1) = sum(C(:))
end

% Closed form formula
counts = N_values.*2.^(N_values - 1);
