% logarithmic distribution
%
% mikael.mieskolainen@cern.ch, 2019

function value = logpdf(k,p)

value = -1./(log(1-p)) .* (p.^k ./ k);

end