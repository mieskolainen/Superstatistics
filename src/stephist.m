% Histogram plot
%
% Plot histogram like stairs() plot but with correct bin centering, which
% stairs() function of Matlab does not have. Stairs() shifts bins
% by one, because it is meant to plot piece-wise functions,
% not histograms per se.
%
% INPUT:    x = Bin centers vector
%           n = Bin values
%
% mikael.mieskolainen@cern.ch, 2019

function h = stephist(x, n, varargin)

x = x(:)';
n = n(:)';

h = stairs([x(1)-(x(2)-x(1))/2 x-(x(2)-x(1))/2 x(length(x))+(x(2)-x(1))/2],[0 n 0], varargin{:});

end
