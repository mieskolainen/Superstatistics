% Shannon entropy
%
% mikael.mieskolainen@cern.ch, 2019

function S = shannonentropy(p)

S = -sum(p(p > 0).*log2(p( p > 0)));

end