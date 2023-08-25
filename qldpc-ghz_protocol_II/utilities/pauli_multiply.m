function [Eout] = pauli_multiply(E1, E2)
% Function to find the product of two Pauli matrices
% E1 = s1 * E(a,b) and E2 = s2 * E(c,d);
% see, for example, https://arxiv.org/abs/1907.00310
% for the E notation

% The first n bits of E1 (resp. E2) give a (resp. c),
% the next n bits give b (resp. d), and
% the last index gives the complex scalar s1 (resp. s2)

% Author: Narayanan Rengaswamy (June 28, 2021)

n = (length(E1)-1)/2;

s1 = E1(end);  % scalar in {1,-1,i,-i} for first Pauli E1
s2 = E2(end);  % scalar in {1,-1,i,-i} for second Pauli E2
sout = s1 * s2;

a = E1(1,1:n);
b = E1(1,(n+1):2*n);  % E1 = E(a,b)
c = E2(1,1:n);
d = E2(1,(n+1):2*n);  % E2 = E(c,d)

sout = sout * 1i^(b*c' - a*d') * (-1)^((a+c) * (b.*d)' + (a.*c) * (b+d)');

Eout = [mod(a+c,2), mod(b+d,2), sout];

end