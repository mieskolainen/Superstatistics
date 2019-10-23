% Formula by Holst, 1980
% http://blog.thegrandlocus.com/static/misc/Holst_1980.pdf
% 
% n stick length x in magnitude: S(1) < S(2) < ... < S(j) < ... < S(n)
%
% P(S_{j} <= x)
%
% mikael.mieskolainen@cern.ch, 2019

function P = sticklength(x, n)

P = zeros(n, length(x));

for j = 1:n
    
for mu = 0:j
    in = 0;
    for nu = 0:n-mu
        in = in + (-1)^nu * nchoosek(n - mu, nu) .* max(0, (1 - (mu + nu)*x)).^(n-1);
    end
    P(j,:) = P(j,:) + nchoosek(n, mu) * in;
end

end

end