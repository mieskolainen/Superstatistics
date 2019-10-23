% Matrix A (inclusion-exclusion matrix) for the compound Poisson inversion
%
% Default bit order follows the order of createCBM() function
%
% input   :       N  =  Vector space dimension
%            reverse =  Reverse the bit order (default false)
% output  :       A  =  Output matrix of size 2^N -1 x 2^N -1
%                       (null/zero vector excluded))
%
% mikael.mieskolainen@cern.ch, 2019

function A = amat(N, reverse)

if (nargin == 1)
    reverse = false;
end

A = zeros(2^N-1);

C = createCBM(N, reverse);
C = C(2:end,:);

% All combinations
ind = 1;
for k = 1:N
   
   M = nchoosek(1:N, k);
   
   for j = 1:size(M,1)
       indices = M(j,:);
       c = C(:,indices);
       
       % Create OR between vectors
       output = zeros(2^N-1, 1);
       for m = 1:size(c,2)
          output = output | c(:,m); 
       end
       
       A(ind,:) = output;
       ind = ind + 1;
   end
end

end