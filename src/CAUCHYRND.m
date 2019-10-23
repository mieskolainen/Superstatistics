% Cauchy random variables
%
% mikael.mieskolainen@cern.ch, 2019

function val = CAUCHYRND(x0,gamma,n)

val = x0 + gamma*tan(pi*(rand(n,1) - 1/2));

end
