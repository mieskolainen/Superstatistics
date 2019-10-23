% Non-linear phase portraits
%
% in (x,dx/dt) phase space
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;

% Time
t = linspace(0, 10, 1e3);

% Set of initial conditions [x,dx/dt]
X0    = [3,   7;
         1,   0;
         -2, -2];

% Random uniform initial conditions
N     = 15;
scale = 6;
X0    = -scale + (scale - -scale)*rand(N,2);

% Loop over initial conditions
% plot x,y lines
plot(linspace(-scale, scale, 4), zeros(4,1), 'k-', 'linewidth', 0.5); hold on;
plot(zeros(4,1), linspace(-scale, scale, 4), 'k-', 'linewidth', 0.5);

for k = 1:size(X0,1)
    % Solve the system numerically for initial conditions
    [T{k},X{k}] = ode45( @F, t, X0(k,:) );
end

% Plot out trajectories
for k = 1:size(X0,1)
    plot(X{k}(:,1), X{k}(:,2)); hold on;
end

% Plot starting point
ax = gca; ax.ColorOrderIndex = 1; % Restart color order
for k = 1:size(X0,1)
    plot(X{k}(1,1), X{k}(1,2), '.', 'Markersize', 7);
end

% This is broken somehow
%u = diff(X(:,1));
%v = diff(X(:,2));
%quiver(X(1:3:end-1,1), X(1:3:end-1,2), u(1:3:end), v(1:3:end));

axis square;
xlabel('$x$','interpreter','latex');
ylabel('$dx/dt$','interpreter','latex');
axis([-scale scale -scale scale]);

