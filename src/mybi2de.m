% Binary vector to decimal
%
% right-msb
% (0,0) = 0, (1,0) = 1, (0,1) = 2, (1,1) = 3
%
% left-msb
% (0,0) = 0, (0,1) = 1, (1,0) = 2, (1,1) = 3
%
% ------------------------------------------------------------------------
%
% Input:  X = N x M  (N = number of bits x M number of samples)
%        order = 'left-msb' or 'right-msb'
%
% Output: Y = M x 1  ( decimal numbers)
%
% mikael.mieskolainen@cern.ch, 2019

function Y = mybi2de(X, order)

    N = size(X,1);
    
    % Output in these
    Y = zeros(size(X,2),1);
    ii = (1:N)';

    if (strcmp(order,'right-msb'))
        for line = 1:size(X,2)
            %vec = X(:,line);
            % Binary expansion
            % result = 0;
            %for i = 1:N
            %    result = result + vec(i)*2^(i-1);
            %end
            %Y(line) = result;
            Y(line) = sum(X(:,line).*2.^(ii-1));
        end
    elseif (strcmp(order,'left-msb'))

        for line = 1:size(X,2)
            vec = X(:,line);
            % Binary expansion
            % result = 0;
            %for i = 1:N
            %    result = result + vec(i)*2^(d-i);
            %end
            %Y(line) = result;
            Y(line) = sum(vec.*2.^(N-ii));
        end
    end
end
