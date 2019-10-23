% Multivariate Hypergeometric Distribution
%
% mikael.mieskolainen@cern.ch, 2019

function p = mhypergeom(K,k)

N = sum(K);
n = sum(k);

% prod_{i=1}^C (K_i k_i) / (N n)
%
% C is the number of different options

val = 1;
for i = 1:C
   val = val * nchoosek(K(i), k(i)); 
end

p = val / nchoosek(N,n);

end