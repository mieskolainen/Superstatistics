% Test symbolic representations
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

N = 2;

% Construct the symbolic representations
[p,m,s,X,EX] = symbolicrep(N);


%% Test Kronecker representations

Nevents = 1e4; N = 2; p = [0.1 0.3 0.6]';
ZETA = zetamat(N);
INVZETA = inv(ZETA);
m = ZETA*[0; p];

% Sample events according to p
S = randsample(1:2^N-1, Nevents, true, p);

% Create binary event vector sample
C = createCBM(N); C = C(2:end,:)';
E = zeros(N,Nevents);
for i = 1:Nevents
   E(:,i) = C(:,S(i));
end

% Calculate Kronecker representations on a finite data sample
pvec = zeros(2^N,1);
mvec = zeros(2^N,1);
svec = zeros(2^N,1);
for i = 1:Nevents
    X1 = [1; E(1,i)];
    X2 = [1; E(2,i)];
    pvec = pvec + INVZETA*kron(X2,X1);
    mvec = mvec + kron(X2,X1);
    
    Y1 = [1; E(1,i) - mean(E(1,:))];
    Y2 = [1; E(2,i) - mean(E(2,:))];
    svec = svec + kron(Y2,Y1);
end
pvec = pvec / Nevents; % Expectation value
mvec = mvec / Nevents;
svec = svec / Nevents;

fprintf('p = [Truth Sampled] \n');
[[0; p], pvec]
fprintf('m = [Truth Sampled] \n');
[m, mvec]
fprintf('svec \n');
svec

% Calculate covariance matrix of the sample
cov(E')


%% Test spectral decomposition of matrix LAMBDA
close all;

%LAMBDA = sym(amat(5));
LAMBDA = amat(10);

[eigvecs,eigvals] = eig(LAMBDA);

plot(real(diag(eigvals)), imag(diag(eigvals)), '.')

