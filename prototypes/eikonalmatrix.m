% Eikonal matrix eigenvalues
%
% For references,
% see e.g http://inspirehep.net/record/1351489/files/Thesis-2014-Roehr.pdf
%         https://arxiv.org/pdf/hep-ph/0007359.pdf
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

% Coupling (real)
syms g z x y

mode = 2;

% Coupling matrix I (no SD with this at all)
if (mode == 1)
g = [z, x - 1i*y;
     x + 1i*y, -z];
end

% Coupling matrix II
if (mode == 2)
g = [z, x - 1i*y;
     x + 1i*y, z];
end

% Coupling matrix III
if (mode == 3)
g = [z, x;
     x, z];
end

% Coupling matrix IV
if (mode == 4)
g = [1, x;
     x, 1];
end

% Eigenmatrix
[eigvec,eigval] = eig(g)

% Dynamical functions
TR  = trace(g)
DET = det(g)

% Kronecker products
% g(X)g
 
gg = kron(g,g)
[S0,D0] = eig(gg);
diag(D0)


%% Triple Pomeron
% g^2(X)g

g2g = kron(g^2,g)
[S1,D1] = eig(g2g);
diag(D1)


%% Triple Pomeron

gg2 = kron(g,g^2)
[eigvec,eigval] = eig(gg2)


%% Double Pomeron

g2g2 = kron(g^2,g^2)
[eigvec,eigval] = eig(g2g2)


%% Check pairwise commutations
% (if commuting -> can be simultaneously diagonalized)

% Commutator
com = @(a,b) a*b - b*a;

com(gg,   g2g)
com(gg,   gg2)
com(gg,  g2g2)
com(g2g, g2g2)
com(gg2, g2g2)

