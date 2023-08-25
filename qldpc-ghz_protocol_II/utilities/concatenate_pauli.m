function [E] = concatenate_pauli(E1, E2)
% Function to concatenate the binary representations of two Pauli matrices
% E1 and E2 acting on different sets of qubits

% Each n-qubit Pauli is a vector of length (2n+1), where the first n
% indices give the X-component, next n indices give the Z-component, and
% the last index contains the overall phase (+1/-1/1i/-1i)

% Without the phase, the vector is generally interpreted as a Hermitian
% Pauli, although this does not affect the function here

% Author: Narayanan Rengaswamy (July 9, 2021)

n1 = (length(E1)-1)/2;
n2 = (length(E2)-1)/2;

E = [ E1(1,1:n1) , E2(1,1:n2) ,  E1(1,(n1+1):2*n1) , E2(1,(n2+1):2*n2) ,  E1(1,end)*E2(1,end) ];

end