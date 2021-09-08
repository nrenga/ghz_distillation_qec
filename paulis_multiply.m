function [Eout] = paulis_multiply(Ein)
% Function to multiply many Pauli matrices which are stored as rows of the
% matrix Ein: each row has (2n+1) entries, the first n giving the
% X-component of the Pauli, the next n giving the Z-component, and the last
% entry giving the overall phase (+1/-1/1i/-1i)

% Author: Narayanan Rengaswamy (July 7, 2021)

n = (size(Ein,2) - 1)/2;
Eout = [ zeros(1,2*n), 1 ];
for i = 1:size(Ein,1)
    Eout = pauli_multiply(Eout, Ein(i,:));
end

end