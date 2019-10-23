% Pileup inversion nonlinear cost function
%
% Any number of unknowns, up to arbitrary (but truncated) Poisson order
%
% Input:    x  = unknown rate vector (2^d-1 x 1)
%
% Output  chi2 = chi2 value
%            Z =  
%            W = 
% mikael.mieskolainen@cern.ch, 2019

function [chi2,Z,W] = nonlinear(x)

% Parameters
global loops;
global nu;
global C;
global Nevents;

% Data
global y;

% Maximum Poisson order
order = length(loops);

% Number of detectors, length of x is 2^N-1
N = log2(length(x)+1);

% Poisson probabilities
P = poisspdf(1:order, nu);
P = P(:) / sum(P);

% Poisson corrections matrix
W = zeros(2^N-1, order);

% Poisson orders k = 1,2,3,...
for k = 1:order
    
    % Fiducial combination weights per order
    wc = zeros(2^N-1, 1);
    
    % Get pile-up combinatorics matrix for this Poisson order
    V = loops{k};
    
    list = cell(2^N-1,1);
    
    % Go through all pile-up combinations, and assign the target
    for v = 1:size(V,1)    
        
        % Gather all overlapping fiducials of this
        comb = zeros(1,N);
        for j = 1:k
            comb = comb + C(V(v,j),:);
        end
        
        % Make it binary
        comb = comb & ones(size(comb));
        
        % Now turn this into fiducial ID. Match is unique because
        % every combination is unique (mutually exclusive)
        dec = mybi2de(comb, 'right-msb');
        list{dec}{end+1} = V(v,:);
    end

    % Multinomial distribution probabilities
    % Create weights
    for i = 1:2^N-1
        %fprintf('%d & ', length(list{i}));
        for j = 1:length(list{i})
            
            % Picking combination
            counts = list{i}{j};
            
            % Count multiplicities
            X = zeros(2^N-1,1);
            for z = 1:2^N-1
               X(z) = sum(counts == z); 
            end
            
            % Multinomial probability (numerical accuracy might be problem
            % here due to product of small numbers -> double underflow)
            %Pm = factorial(k) / (prod(factorial(X))) * prod( x .^ X );
            Pm = factorial(k) / (prod(factorial(X))) * exp(sum(log( x .^ X )));
            
            % Add the probability
            wc(i) = wc(i) + Pm;
        end
        %fprintf('\n');
    end
    %fprintf('\n');

    % Add this column
    W(:,k) = wc;
end

% Finally generate the observation
Z = W * P;

% Chi^2 Cost with event counts and sqrt(N) errors
chi2  = sum((y - Z).^2);
%chi2 = sum((y*Nevents - Z*Nevents).^2 ./ (y*Nevents) );

end

%% FIXED VERSION FOR REFERENCE

%{
function [chi2,Z] = analytic(x)

% Parameter
global nu;
global N_events;

% Data
global y;

% Poisson probabilities
P = [poisspdf(1,nu) poisspdf(2,nu) poisspdf(3,nu) poisspdf(4,nu) poisspdf(5,nu)]; P = P / sum(P);

% Joint prob matrix for k = 2
J = zeros(3,3);
i = 0;
for a = 1:3
    for b = 1:3
        if (b >= a)
            J(a,b) = x(a)*x(b);
        end
    end
end

J = J / sum(J(:));

% Joint prob matrix for k = 3
S = zeros(3,3,3);
for a = 1:3
    for b = 1:3
        for c = 1:3
            if ((c >= b) && (b >= a))
               S(a,b,c) = x(a)*x(b)*x(c);
            end
        end
    end
end
S = S / sum(S(:));

% Joint prob matrix for k = 4
D = zeros(3,3,3,3);
for a = 1:3
    for b = 1:3
        for c = 1:3
            for d = 1:3
                if ((d >= c) &&(c >= b) && (b >= a))
                   D(a,b,c,d) = x(a)*x(b)*x(c)*x(d);
                   i = i + 1;
                   fprintf('D(%d,%d,%d,%d)\n',a,b,c,d);
                end
            end
        end
    end
end
D = D / sum(D(:));

% Joint prob matrix for k = 5
R = zeros(3,3,3,3,3);
for a = 1:3
    for b = 1:3
        for c = 1:3
            for d = 1:3
                for e = 1:3
                    if ((e >= d) && (d >= c) &&(c >= b) && (b >= a))
                       R(a,b,c,d,e) = x(a)*x(b)*x(c)*x(d)*x(e);
                       i = i + 1;
                       fprintf('R(%d,%d,%d,%d,%d)\n',a,b,c,d,e);
                    end
                end
            end
        end
    end
end
R = R / sum(R(:));


% One interaction:     k = 1
Z = [P(1)*x(1);
     P(1)*x(2);
     P(1)*x(3)];

% Two interactions:    k = 2
Z = Z + ...
    [P(2)*J(1,1);
     P(2)*J(2,2);
     P(2)*(J(3,3) + J(1,2) + J(1,3) + J(2,3))];

% Three interactions:  k = 3
Z = Z + ...
    [P(3)*S(1,1,1);
     P(3)*S(2,2,2)
     P(3)*(S(3,3,3) + S(1,1,2) + S(1,1,3) + S(1,2,2) + S(1,2,3) + S(1,3,3) + S(2,2,3) + S(2,3,3)) ];

% Four interactions:   k = 4
 Z = Z + ...
     [P(4)*D(1,1,1,1);
      P(4)*D(2,2,2,2);
      P(4)*(D(3,3,3,3) + D(1,1,1,1) + D(1,1,1,3) + D(1,1,2,2) + D(1,1,2,3) + D(1,1,3,3) ...
          + D(1,2,2,2) + D(1,2,2,3) + D(1,2,3,3) + D(1,3,3,3) + D(2,2,2,3) + D(2,2,3,3) + D(2,3,3,3)) ];

% Five interactions:   k = 5
Z = Z + ...
     [P(5)*R(1,1,1,1,1);
      P(5)*R(2,2,2,2,2);
      P(5)*(R(3,3,3,3,3) + R(1,1,1,1,2) + R(1,1,1,1,3) + R(1,1,1,2,2) + R(1,1,1,2,3) + R(1,1,1,3,3) + R(1,1,2,2,2) +  ...
            R(1,1,2,2,3) + R(1,1,2,3,3) + R(1,1,3,3,3) + R(1,2,2,2,2) + R(1,2,2,2,3) + R(1,2,2,3,3) + R(1,2,3,3,3) +  ...
            R(1,3,3,3,3) + R(2,2,2,2,3) + R(2,2,2,3,3) + R(2,2,3,3,3) + R(2,3,3,3,3)) ];

% Chi^2 Cost with event counts and sqrt(N) errors
chi2 = sum((y*N_events - Z*N_events).^2 ./ (y*N_events) )

end
%}