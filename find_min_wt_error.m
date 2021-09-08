function [err_est] = find_min_wt_error(S, syndrome)
% Find the minimum weight Pauli error matching the syndrome
% for the stabilizer code defined by the stabilizer matrix S

% S is a [m x (2n+1)] matrix where the last column is sign (+/- 1);
% syndrome is a binary vector with 1s indicating the rows of
% S with which the actual error anticommutes

% If an error does not match the syndrome, the function automatically eliminates
% all patterns in the stabilizer coset determined by that error

% Author: Narayanan Rengaswamy (July 9, 2021)

n = (size(S,2) - 1)/2;
m = size(S,1);
syndrome = syndrome(:);  % make it a column vector
eliminated = [];

if (all(syndrome == 0))
    err_est = zeros(1,2*n);
    return;
end

% Pauli:  I, X, Y, Z
pauli = [ 0, 1, 1, 0 ; ...
          0, 0, 1, 1 ];

for i = 1:n   % determines weight of the Pauli operator
    qubits = nchoosek(1:n,i);
    for combo = 1:size(qubits,1)
        for j = 1:(4^i - 1)
            paulis = de2bi(j,i,4) + 1;
            err = zeros(2,n);
            err(:,qubits(combo,:)) = pauli(:,paulis);
            err = reshape(err', 1, 2*n);
            if (any(eliminated == bi2de(err)))
                continue;
            else
                if (all(mod(S(:,1:2*n) * fftshift(err,2)', 2) - syndrome == 0))
                    err_est = err;
                    return;
                else
                    eliminated = [eliminated; bi2de(mod(err + mod(de2bi(0:(2^m-1),m,2) * S(:,1:2*n),2), 2))];
                end
            end
        end
    end
end

end