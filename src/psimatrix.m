% Create psi-matrix
% Note that this matrix is rank deficient:
% 
% d  2^d-1  rank(psi)
% ---------------
% 2  3      3
% 3  7      6
% 4  15     10
% 5  31     15
% 6  63     21
% 7  127    28
% 8  255    36
%
% input:    d = dimension
%         
% output:  XI = X'*X
%         PSI = XI.^2 matrix (2^d-1 x 2^d-1)
%
% mikael.mieskolainen@cern.ch, 2019

function [PSI,XI] = psimatrix(d)

X = createCBM(d)';
X = X(:,2:end);     % Remove the first row

% Normalize each vector
for i = 1:size(X,2)
    X(:,i) = X(:,i) / norm(X(:,i));
end

% Create psi-matrix
XI = (X')*X;
PSI = XI.^2;

end