% Condition number of inclusion-exclusion (Mobius inversion) matrices
%
% mikael.mieskolainen@cern.ch, 2019

N_values = 2:10;
c = zeros(length(N_values),1);

for i = 1:length(N_values)
    A = amat(N_values(i));
    c(i) = cond(A);
end

figure;
semilogy(N_values, c);
xlabel('$N$', 'interpreter', 'latex');
ylabel('cond(A)', 'interpreter', 'latex');
set(gca,'fontsize',13);

