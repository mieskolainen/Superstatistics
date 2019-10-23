% Random numbers from logpdf
%
% mikael.mieskolainen@cern.ch, 2019

function values = logrnd(p, n)

MAX = 10; % Maximum cutoff for speed
pdf = logpdf(1:MAX, p);

values = zeros(n,1);
bin = 0;
for i = 1:n
    
    % Pick bin
    while (true)
        bin = randi(MAX);
        if (rand(1) < pdf(bin))
            break;
        end
    end
    values(i) = bin;
end

end