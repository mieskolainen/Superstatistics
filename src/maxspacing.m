% Maximum Spacing statistics
%
% https://link.springer.com/content/pdf/10.1007%2Fs00362-008-0134-3.pdf
%
% mikael.mieskolainen@cern.ch, 2019

function P = maxspacing(x,n)

P = zeros(n+1,1);

for m = 1:n+1
    summa = 0;
    for i = m:n+1
        summa = summa + (-1)^(i-1)/(n+2-i) * nchoosek(n,i-1) * (1-x*(n+2-i))^n;
    end
    summa = summa * (-1)^n * (n+1);
    
    P(m) = summa;
end

end