% Create M matrix
% Calculate Zeta upper triangular matrix for Galois field GF(2^N)
%
% [1 -1
% [0 1] to Kronecker product ^ N
%
% see also zetamat()
%
% mikael.mieskolainen@cern.ch, 2019

function M = mmat(N)

L = [1 -1;
     0  1];

M = L;
for i = 1:N-1
    M = kron(M,L);
end

end