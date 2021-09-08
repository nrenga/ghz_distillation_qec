function [Sout, first_anticommute] = stabilizer_formalism_msmt(Sin, E)
% Function that updates the stabilizer generators in Sin after the
% measurement of the Pauli operator E, based on the stabilizer formalism
% for measurements (https://arxiv.org/abs/quant-ph/9807006)

% The first 2n bits of E give the binary representation of the Pauli matrix
% and the last index provides the measurement result, i.e., sign, 1 or -1;
% each row of Sin is a stabilizer generator for the input state and has the 
% same format as E

% Author: Narayanan Rengaswamy (June 29, 2021)

n = (length(E)-1)/2;
first_anticommute = 0;
Sout = Sin;

for i = 1:size(Sout,1)
    if (symp_inn_pdt(Sout(i,1:2*n), E(1,1:2*n)) == 1)
        if (first_anticommute == 0)
            first_anticommute = i;
            Sout(i,:) = E;   % last index gives syndrome (sign)
        else
            Sout(i,:) = pauli_multiply(Sin(i,:), Sin(first_anticommute,:));
        end
    end
end

end