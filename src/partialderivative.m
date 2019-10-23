% d\mathbf{y}/du partial derivatives
% 
% mikael.mieskolainen@cern.ch, 2019
function dydu = partialderivative(N)

A = sym(amat(N));
p = sym('p',[2^N-1 1]);
u = sym('u');

dydu = inv(A) * ( ( exp(-u*A*p)-1) / (exp(u)*(exp(-u)-1)^2) - (A*p) .* exp(-u*A*p)/(exp(-u)-1) );

% Tested with:
% dtest = diff(inv(A) *( exp(-u*A*p) - 1) / (exp(-u)-1), u)
% simplify(dydu - dtest)

end