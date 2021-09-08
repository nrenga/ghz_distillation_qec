function [S] = find_binary_symmetric_matrix(A,B)
% Function to find a binary symmetric matrix S satisfying A*S = B (mod 2)
% Vectorize: A*S*I = B  => kron(I^T,A) * S(:) = B(:) subject to
%                S = S' => S(:) = Perm * S(:),
% where Perm is the permutation matrix that converts S(:) to (S')(:)

% vec(A*B*C) = (C^T \otimes A) * vec(B)

% Author: Narayanan Rengaswamy (July 2, 2021)

n = size(A,2);
S = reshape((1:n^2)',n,n);
I = eye(n^2);
Perm = I(reshape(S',1,n^2),:);

Aext = [ kron(eye(n),A) ; mod(I - Perm,2) ];
b = [ B(:) ; zeros(n^2,1) ] ;
s = gflineq(Aext,b,2);

S = reshape(s,n,n);

end