% Create Mobius inversion matrices in C++ format
%
% mikael.mieskolainen@cern.ch
clear; close all;
addpath src

NLIST = 1:8;

fp = fopen('matrix.cc','w');

for k = 1:2 % Matrix and its inverse
    
    for N = NLIST % Loop over all N
        
        LAMBDA = amat(N, true);
        if (k == 1)
            fprintf(fp, 'std::vector<std::vector<double> > A%d = \n{', N);
        else
            fprintf(fp, 'std::vector<std::vector<double> > invA%d = \n{', N);
            LAMBDA = inv(LAMBDA);
        end
        
        for i = 1:size(LAMBDA,1)
            fprintf(fp, '{');
            for j = 1:size(LAMBDA,2)
                if (j < size(LAMBDA,2))
                fprintf(fp, '%d, ', LAMBDA(i,j));
                else
                    fprintf(fp, '%d', LAMBDA(i,j));
                end
            end
            if (i < size(LAMBDA,1))
                fprintf(fp, '},\n');
            else
                fprintf(fp, '}};\n');
            end
        end
    end
end

fclose(fp);
