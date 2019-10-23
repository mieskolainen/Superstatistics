% AGK cutting rules toy simulation
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

N = 1e5;

% Free parameters
tot   = 107;
DELTA = 0.10;

s = 13000^2;
lambda = log(s^DELTA);
mMAX = ceil(2*lambda);

wmat = zeros(mMAX, mMAX + 1);
amat = zeros(mMAX, mMAX + 1);
smat = zeros(mMAX, mMAX + 1);

for i = 1:N

    % Number of exchanges according to Poisson
    while(true)
        m = poissrnd(lambda);
        if (m > 0 && m <= mMAX)
            break;
        end
    end
    
    % Number of cuts flat [0,1,..,m]
    k = randi(m+1)-1;
    
    % weight and sign
    w = AGK(m,k);
    
    wmat(m,k+1) = wmat(m,k+1) + w;
    amat(m,k+1) = amat(m,k+1) + abs(w);
    smat(m,k+1) = smat(m,k+1) + sign(w);
end

wmat = wmat / N
amat = amat / N;
smat = smat / N;

X = [];
for u = 1:mMAX
    for nu = 0:u
        X(u,nu+1) = AGK(u,nu);
    end
end

sum(wmat(:))
E = amat .* sign(X)

sum(E(:))

%% Binomial vs Poisson comparison
close all;

nbdpdf = @(n, avgN, k) gamma(n+k) ./ (gamma(n+1) .* gamma(k)).*(avgN ./ (k + avgN)).^n .* (k / (k + avgN)).^k;

avgN = 2*log(300^2)
k = avgN*4;

x = 1:3*avgN;

%val = nbdpdf(n, avgN, k);

p = 0.3;
N = round(avgN/p);
val = binopdf(x, N, p)

plot(x, val); hold on;
plot(x, poisspdf(x, avgN), 'r');

legend('Binomial','Poisson');

%set(gca,'yscale','log');
