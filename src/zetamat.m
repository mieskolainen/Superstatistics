% Calculate Zeta upper triangular matrix for Galois field GF(2^N)
% Note that the inverse of Zeta is the so-called Mobius matrix.
%
% Input:       N  =  dimension of the vector space
% Output:    zeta =  Zeta matrix (2^N x 2^N)
%
% The zeta matrix for a vector space V is
%       zeta(S,T) = 1, if S \subset T
%                   0, otherwise
%
% Zeta is also here the Sierpinski triangle matrix
%
% And can be calculate with Kronecker tensor product
% 
% X = [1 -1; 
%      0  1];
% kron(X,X)  ( same as inv(zetamat(2)) )
% 
% mikael.mieskolainen@cern.ch, 2019

function zeta = zetamat(N)

C = createCBM(N);
zeta = zeros(2^N,2^N);

for i = 1:size(zeta,1)
    for j = i:size(zeta,2)
        if (sum( (C(i,:) | C(j,:)) == C(j,:)) == N )
        zeta(i,j) = 1;
        end
    end
end

end