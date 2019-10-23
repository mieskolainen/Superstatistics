% Symbolic derivative test
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src


%%
syms u x

fun = int(taylor((exp(-u) - exp(-(2*u)/3))*exp(-1i*x*u), u, 'ExpansionPoint', 10, 'Order', 3), u);

xval = linspace(1,10,1e2);
vals = zeros(size(xval));

for i = 1:length(xval)
    vals(i) = double( subs(fun, [x; u], [xval(i); 10]) );
end

plot(xval, real(vals))
    
%% Loop over
for N = 3

A = sym(amat(N));
p = sym('p',[2^N-1 1]);
u = sym('u');

% Map p -> y
y = inv(A) *( exp(-u*A*p) - 1) / (exp(-u)-1);

% Shannon entropy
S = simplify(-sum(y .* log(y)));

% C.) Create basis vector input
C = createCBM(N);
C = C(2:end,:);
CSUM = sum(C,2);
ind = find(CSUM == 1);
input = zeros(2^N-1,1);
input(ind) = 1;
input = input / sum(input);
pval = input(:)';

% Create symbolic variables
for i = 1:2^N-1
    eval(sprintf('p%d = sym(''p%d'');', i, i)); 
end

% mu-values
uval = linspace(-100,30,1e2);
vals = zeros(2^N-1,length(uval));

% Create expression
str = 'val = subs(y, [';
for i = 1:2^N-1
    str = [str, sprintf('p%d ', i)];
end
str = [str, '], [pval])'];

eval(str)

end