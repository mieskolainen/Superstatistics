% Random numbers from Poisson distribution
%
% mikael.mieskolainen@cern.ch, 2019

function values = POISRND(lambda, n)

MAX = ceil(lambda * 3); % Maximum cutoff for speed
pdf = poisspdf(0:MAX-1, lambda);

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
    values(i) = bin-1;
end

end