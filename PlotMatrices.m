% Plot combinatorial matrices
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src

N = 7;
LAMBDA = amat(N);
LAMBDAINV = inv(LAMBDA);

f1 = figure;
imagesc(LAMBDA); axis square; colormap(hot)

    filename = sprintf('../figs/A7x7.pdf');
    print(f1, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

f2 = figure;
imagesc(LAMBDAINV); axis square; colormap(hot)

    filename = sprintf('../figs/Ainv7x7.pdf');
    print(f2, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));


%%

N = 7;
ZETA = zetamat(N);
ZETAINV = inv(ZETA);

f3 = figure;
imagesc(ZETA); axis square; colormap(hot)

    filename = sprintf('../figs/ZETA%dx%d.pdf', N, N); 
    print(f3, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

f4 = figure;
imagesc(ZETAINV); axis square; colormap(hot)

    filename = sprintf('../figs/ZETAinv%dx%d.pdf', N, N); 
    print(f4, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));


%% "Fish skeleton matrix"

f5 = figure;
LL = LAMBDA*LAMBDA';
imagesc(LL); axis square;

    filename = sprintf('../figs/lambdalambdaT.pdf'); 
    print(f5, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

f6 = figure;
ZZ = ZETA*ZETA';
imagesc(ZZ); axis square; colormap(hot);

    filename = sprintf('../figs/zetazetaT.pdf'); 
    print(f6, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));


%%

f7 = figure;
imagesc(LAMBDAINV*LAMBDAINV'); axis square; colormap(hot)
pdfcrop(0.69, 0.91);

    filename = sprintf('../figs/lambdainvlambdainvT.pdf'); 
    print(f7, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

f8 = figure;
imagesc(ZETAINV'*ZETAINV); axis square; colormap(hot)
pdfcrop(0.69, 0.91);

    filename = sprintf('../figs/zetainvTzetainv.pdf'); 
    print(f8, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));


%%

f9 = figure;
zeta = zetamat(N);
%zeta = zeta(2:end, 2:end);

subplot(2,4,1); imagesc(zeta); colormap hot; axis square; title('$\zeta$','interpreter','latex');
subplot(2,4,2); imagesc(inv(zeta)); axis square; title('$\zeta^{-1}$','interpreter','latex');

A = amatfull(N);
A_inv = inv(A);
subplot(2,4,3); imagesc(A); axis square; title('$A$','interpreter','latex');
subplot(2,4,4); imagesc(A_inv); axis square; title('$A^{-1}$','interpreter','latex');

M = (A*inv(zeta))';
M_inv = ~fliplr(inv(double(M)));
subplot(2,4,5); imagesc(M); axis square; title('$M = (\zeta A^{-1})^T$','interpreter','latex');
subplot(2,4,6); imagesc(M_inv); axis square; title('$M^{-1}$','interpreter','latex');

%colormap gray;
% M*zeta = A
% M = A*inv(zeta)

% zeta = M*A
% M = zeta * inv(A)

    
%% Adjacency matrix
close all;

N = 6;
A = hypercube(N);

f10 = figure;
imagesc(0:2^N-1, 0:2^N-1, A); colormap(hot);
xticks([0:8:2^N]);
yticks([0:8:2^N]);
axis square;
    
    filename = sprintf('../figs/adjmatrixN%d.pdf', N); 
    print(f10, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

% Degree matrix (each vertex is connected by N edges)
D = eye(2^N)*N; % Same as diag(sum(A))

% Laplacian matrix L = D - A
L = D - A;

% Symmetric Laplacian matrix
Lsym = eye(2^N) - D^(-1/2) * A * D^(-1/2);

imagesc(Lsym);


%%
% Resistance distance
% http://mathworld.wolfram.com/ResistanceDistance.html
Gamma = L + 1/(2^N);
Omega = zeros(2^N);
invGamma = inv(Gamma);
for i = 1:2^N
    for j = 1:2^N
        Omega(i,j) = invGamma(i,i) + invGamma(i,j) - 2*invGamma(i,j);
    end
end

figure;
imagesc(Omega);
axis square;

% Graph connections
GR = graph(A);

fg = figure;
plot(GR,'layout','circle','NodeLabel',0:2^N-1);
axis square;
axis off;


    filename = sprintf('../figs/graphN%d.pdf', N); 
    print(fg, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', filename, filename));

