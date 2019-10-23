% n Powerlaw random numbers from [x0,x1], alpha the exponent
%
% mikael.mieskolainen@cern.ch, 2019

function val = POWERF(x0,x1,alpha,n)

val = (( (ones(n,1)*x1).^(alpha+1) - (ones(n,1)*x0).^(alpha+1)).*rand(n,1) ...
       + (ones(n,1)*x0).^(alpha+1)).^(1./(alpha+1));

end