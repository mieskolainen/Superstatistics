% Determinants, traces, sums etc.
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

dets   = [];
traces = [];
sums   = [];
permanents = [];

N_values = 2:8;
mode = 0;  % Forward or inverse

for N = N_values
    if (mode == 0)
        LAMBDA = amat(N);
    else
        LAMBDA = inv(amat(N));
    end

    if (N < 3)
    (-1)^sym(LAMBDA)
    end
    
    dets(end+1) = det(LAMBDA);
    traces(end+1) = trace(LAMBDA);
    sums(end+1) = sum(LAMBDA(:));
    %permanents(end+1) = permanent_mat(LAMBDA);

    [vv,xx] = eig(LAMBDA);
    dd = diag(xx);
    plot(real(xx), imag(xx), 's');
end

%%
stem(N_values, dets);
xlabel('$N$','interpreter','latex');
ylabel('det($\Lambda$)','interpreter','latex');
