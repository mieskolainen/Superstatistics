% Symbolic Jacobians
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath ../src

%% Create symbolic matrices

N = 2; % Vector space dimension
mu = sym('mu');

reverse = false;
A = amat(N, reverse);
A = sym(A);
A_inv = inv(A);

% Choose direction here (0 or 1)
for direction = 0

if (direction == 0)

fprintf('\n\n ** Forward equation ** \n');

%% Forward equation
p = sym('p',[2^N-1 1]);
y = A_inv * (exp(-mu*A*p) - 1) / (exp(-mu)-1);

if (N >= 2)
%% Create unitarity string p_n = 1 - (p1+p2+...+p_n-1) (Forward)

unitarity = '1 - ';
for i = 1:2^N-3
    unitarity = [unitarity, sprintf('p%d - ', i)];
end
unitarity = [unitarity sprintf('p%d', size(A,1)-1)];
fprintf('Unitarity:\n');
disp(unitarity);


%% Substitute unitarity (Forward)

y = subs(y, sprintf('p%d', size(A,1)), str2sym(unitarity));
%p = simplify(p, 25)
y = simplify(y, 'IgnoreAnalyticConstraints', true, 'Step', 25); % ignore analytic constraint gives: log(exp(-mu)) = -mu

end

y

%% Calculate Jacobian matrix (Forward)

last_index = max(length(y)-1, 1); % Take care of N = 1 case

% Add also differentiation wrt mu here
J = simplify( jacobian(y, [p(1:last_index)] ), 'Step', 25)
rank(J)



else

fprintf('\n\n ** Inverse equation ** \n');


%% Inverse equation
y = sym('y',[2^N-1 1]);
p = A_inv*log((exp(-mu)-1)*A*y+1)/-mu;


if (N >= 2)
    
%% Create unitarity string y_n = 1 - (y1+y2+...+y_n-1) (Inverse)

unitarity = '1 - ';
for i = 1:2^N-3
    unitarity = [unitarity, sprintf('y%d - ', i)];
end
unitarity = [unitarity sprintf('y%d', size(A,1)-1)];
fprintf('Unitarity:\n');
disp(unitarity);


%% Substitute unitarity (Inverse)

p = subs(p, sprintf('y%d', size(A,1)), unitarity);
%p = simplify(p, 25)
p = simplify(p, 'IgnoreAnalyticConstraints', true, 'Step', 25); % ignore analytic constraint gives: log(exp(-mu)) = -mu

end

p
%% Calculate Jacobian matrix (Inverse)

last_index = max(length(y)-1, 1); % Take care of N = 1 case

% Add also differentiation wrt mu here
J = simplify( jacobian(p, [y(1:last_index)] ), 'Step', 25)
rank(J)

end

end % For-loop
