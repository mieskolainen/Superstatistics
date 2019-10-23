% Sample from Tsallis pdf
%
% mikael.mieskolainen@cern.ch, 2019

function samples = tsallisrnd(q,T, maxpt, N)

func = @(mt, T, q) (1+(q-1) * mt/T).^(q/(1-q));

% q = 1.13;
% T = 0.084;

% n = 1/(q-1);
% p0 = T/(q-1);

% Pion mass
mpi = 0.139;

% Discretization
bins = 1e3;

ptval = linspace(0, maxpt, bins);
mtval = sqrt(ptval.^2 + mpi^2);

% Evaluate dsigma/dpt |_y=0 (mid rapidity)
dsdpt = ptval .* mtval .* func(mtval, T, q);
maxval = max(dsdpt);

% Acceptance-Rejection
samples = zeros(N,1);
k = 1;
while (true)
    ID = randi(bins);
    if (rand(1)*maxval < dsdpt(ID))
        samples(k) = ptval(ID);
        k = k + 1;
    end
    if (k == N+1)
        break;
    end
end
end