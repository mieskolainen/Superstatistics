% Matrix A (inclusion-exclusion matrix) for
% combinatorial compound Poisson inversion
%
% Input   :    d  =  Integer dimension
% Output  :    A  =  Output matrix of full size 2^d x 2^d
%
% mikael.mieskolainen@cern.ch, 2019

function A_full = amatfull(d, reverse)

if (nargin == 1)
    reverse = false;
end

A = zeros(2^d-1);

C = createCBM(d, reverse);
C = C(2:end,:);

% All combinations
ind = 1;
for k = 1:d
   
   M = nchoosek(1:d, k);
   
   for j = 1:size(M,1)
       indices = M(j,:);
       c = C(:,indices);
       
       % Create OR between vectors
       output = zeros(2^d-1, 1);
       for m = 1:size(c,2)
          output = output | c(:,m); 
       end
       
       A(ind,:) = output;
       ind = ind + 1;
   end
end

% Add the contribution from 0-vector (basically ones to left most column
% and top most row)
A_full = zeros(size(A,1)+1, size(A,2)+1);
A_full(2:end, 2:end) = A;
A_full(1,1) = 1;

end
