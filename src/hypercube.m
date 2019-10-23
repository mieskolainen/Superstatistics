% Graph theory adjacency matrix of N-hypercube
%
% The code from:
% https://blogs.mathworks.com/cleve/2017/02/20/hypercubes-and-graphs/
%
% mikael.mieskolainen@cern.ch, 2019

function A = hypercube(N)

    pow2 = 2.^(0:N-1);
    f2b = @(j) floor(rem(j./pow2, 2)); % to binary
    b2f = @(b) b*pow2';                % to int
    n   = 2^N;
    A   = zeros(n, n, 'logical');
    
    for j = 0:n-1
        % Column indices to binary
        b = f2b(j);
        
        % Flip bits for the row indices
        for i = 1:N
            b(i) = ~b(i);
            k    = b2f(b);
            b(i) = ~b(i);
            A(k+1,j+1) = 1;
        end
    end
end