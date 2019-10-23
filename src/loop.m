% Any dimensional for-for-for-...-loop generation
%
% Input: n  =  2^N-1, # number of non-zero fiducial rates, e.g. 7
%        K  =  Order, e.g. 5 (such as Poisson truncation order)
%
% mikael.mieskolainen@cern.ch, 2019

function VM = loop(n,K)

V  = ones(1,K);
go = true;

% Count the multinomial multiplicity and make the output matrix
multiplicity = nchoosek(K+n-1,n-1);
VM = zeros(multiplicity, K);

% First one is filled with ones, so start from index 2
j = 1;

while (go)
   
    % Check the save condition
   save = true;
    for i = K:-1:2
        if (V(i) >= V(i-1))
            save = true;
        else
            save = false;
            break;
        end
    end
    
    if (save)
        VM(j,:) = V;
        j = j + 1;
    end
    
    % Increment
    V(K) = V(K) + 1;
    if V(K) > n
        V(K) = 1;
        go   = false;
        for i = K-1:-1:1
            V(i) = V(i) + 1;
            if (V(i) <= n)   % i-th counter not at the limit
                go = true;
                break;       % Leave "for i" loop
            end
            V(i) = 1;        % Reset i-th counter 
        end
    end
end

end