% Print matrix elements in Mathematica (ascii) format
%
% input: A = matrix
%
% mikael.mieskolainen@cern.ch, 2019

function mtomat(A)

fprintf('{');
for i = 1:size(A,1)
    fprintf('{');
    for j = 1:size(A,2)
        if (j < size(A,2))
            fprintf('%d, ', A(i,j));
        else
            fprintf('%d', A(i,j)); 
        end
    end
    if (i < size(A,1))
        fprintf('},');
    else
        fprintf('}');
    end
end
fprintf('} \n');

end