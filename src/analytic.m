% Analytic combinatorial compound Poisson inversion
%
% input:      y  =  Measured normalized rates, vector  (2^N-1 x 1)
%          R_OR  =  "Global" OR rate [0 ... 1]         scalar
%
% output:  p_hat =  Inverse estimate rates    (2^N-1 x 1)
%
% mikael.mieskolainen@cern.ch, 2019

function p_hat = analytic(y, R_OR, LAMBDA, invLAMBDA)

% A and its inverse as global variables for speed
%global LAMBDA;
%global invLAMBDA;

if (nargin == 2)
    % OR Create mapping matrix here
    N = log2(length(y)+1); % Number of observables
    LAMBDA = sparse(amat(N));
    invLAMBDA = inv(LAMBDA);
end

% Inverse
r =  log((exp(-mu(R_OR))- 1)*LAMBDA*y + 1) / -mu(R_OR);
p_hat = invLAMBDA * r;

% The inverse also with polylog(1,z)
%p_hat = LAMBDA \ ( polylog(1,LAMBDA*y*R_OR) / polylog(1,R_OR) );

end

% Rate to Poisson mu-value function
function mu_value = mu(R)
mu_value = -log(1-R);
end
