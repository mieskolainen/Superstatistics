% Symbolic representations
%
% Input:  N = Vector space dimension (>= 2)
%
% Output: p = Natural/fundamental representation (symbolic)
%         m = Moment/correlation representation (symbolic)
%         s = Centralized moment (x - E[x]) representation (symbolic)
%         X = Matrix of symbolic variables [1 x1]' [1 x2]' [1 x3]'
%        EX = Matrix of centralized variables [1 x1-E[x1]]' ... (symbolic)
%
% mikael.mieskolainen@cern.ch, 2019

function [p,m,s,X,EX] = symbolicrep(N)

if (N < 2)
   error('Input needs to be N >= 2');
   return; 
end

% Create zeta and the inverse matrix
ZETA = sym(zetamat(N));
invZ = sym(inv(ZETA));

% Create a symbolic matrix of X variables [1 x1]' [1 x2]' etc.
X = sym(zeros(2,N));
EX = X;
for i = 1:N
   X_this = sym('x',[2 1]); X_this(1) = 1;   X_this(2) = sprintf('x%d',i);
   EX_this = sym('x',[2 1]); EX_this(1) = 0; EX_this(2) = sprintf('Ex%d',i); 
   X(:,i) = X_this;
   EX(:,i) = X_this - EX_this;
end

%X1 = sym('x',[2 1]); X1(1) = 1; X2(2) = 'x1'
%X2 = sym('x',[2 1]); X2(1) = 1; X2(2) = 'x2'
%X3 = sym('x',[2 1]); X3(1) = 1; X3(2) = 'x3'

% Any d recursively
rec = kron(X(:,2),X(:,1));
Erec = kron(EX(:,2),EX(:,1));
if (N >= 3)
    for i = 3:N
        rec = kron(X(:,i),rec);
        Erec = kron(EX(:,i),Erec);
    end
end
% Fixed N = 2 or N = 3
%if (N == 2)
%    p = invZ*kron(X2,X1)
%elseif (N == 3)
%    p = invZ*kron(X3,kron(X2,X1))    
%end

% p and m representations
p = invZ*rec;
m = ZETA*p;
s = Erec;

end