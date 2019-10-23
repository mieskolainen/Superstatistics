% n Uniform random numbers from [a,b]
%
% mikael.mieskolainen@cern.ch, 2019

function val = UNIF(a,b,n)

val =  a + (b - a)*rand(n,1);

end
