% Test AGK (Abramovski-Gribov-Kancheli) rules
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath src

kronsum = @(A,B) kron(A, eye(size(A))) + kron(eye(size(B)), B);

fprintf('\n');
C = zeros(12);

% Bartels-Ryskin matrix
R = [0 0 1/2;
     0 1 0;
     1 0 1/2];

% Number of Pomerons exchanged
for m = 1:size(C,1)
    fprintf('%d ', m);
    % Number of Pomerons cut
    for k = 0:m
        C(m,k+1) = AGK(m,k);
    end
    fprintf('\n');
end

C
C'

for i = 1:size(C,1)
    fprintf('%d & ', i);
    for j = 1:size(C,2)
        fprintf('%d & ', C(i,j)); 
    end
    fprintf('\\\\ \n');
end

%% Plot

% Kronecker product of transition matrices
f1 = figure;
PRO = kron(R, kron(R, kron(R, R)));
imagesc(PRO); axis square; colormap(hot);

filename = sprintf('../figs/AGK_tensorprod.pdf');
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));

% Kronecker sum of transition matrices
f2 = figure;
SUM = kronsum(kronsum(R,R),kronsum(R,R))/4;
imagesc(SUM); axis square; colormap(hot);

filename = sprintf('../figs/AGK_tensorsum.pdf');
print(f2, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


