% Recursive loops
%
% mikael.mieskolainen@cern.ch, 2019

function loop_rec(X,y,n)

if (n >= 1) 
    for i = 1:length(y)
        loop_rec(X, y, n - 1)
    end
else
    for i = 1:size(X,1)
    for j = 1:size(X,2)
       fprintf('%0.5f\n', X([i j]));
    end
    end
end

end