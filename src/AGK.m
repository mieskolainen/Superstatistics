% AGK (Abramovski-Gribov-Kancheli) Cutting Rules Combinatorics
%
% Input:  m = Number of Pomerons exchanged
%         k = Number of Pomerons cut, 0 <= k <= m
%
% -1 : 4 : -2 etc. factors
% 
% Output:  c = AGK factor
%
% mikael.mieskolainen@cern.ch, 2019

function c = AGK(m,k)

c = (-1)^(m-k)*factorial(m)/(factorial(k)*factorial(m-k))*(2^(m-1) - delta(0,k));

end

function out = delta(i,j)

if (i == j)
    out = 1;
else
    out = 0;
end

end