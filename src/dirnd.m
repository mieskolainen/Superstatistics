% Random numbers from the Diricthlet distribution
%
% Example: dirnd(ones(5,1), 1000);
%
% Input:   a  =  Parameter vector controlling the distribution
%          n  =  Number of random vectors
%
% mikael.mieskolainen@cern.ch, 2019

function r = dirnd(a,n)

p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

end