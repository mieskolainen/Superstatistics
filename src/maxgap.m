% Maximum rapidity gap (spacing statistics)
%
% mikael.mieskolainen@cern.ch, 2019

function P = maxgap(nval, beta)

P = zeros(length(nval),1);
for i = 1:length(nval)
    n = nval(i);
    for j = 0:n
        if (j < 1/beta)
            P(i) = P(i) + nchoosek(n+1,j) * (-1)^j * (1 - beta*j)^n;
        end
    end
end

end