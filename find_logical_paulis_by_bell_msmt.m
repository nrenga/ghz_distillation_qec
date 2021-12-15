function [Xbar, Zbar] = find_logical_paulis_by_bell_msmt(Hin)
% Function to compute logical X and Z generators for the stabilizer code
% defined by the parity check matrix Hin, a [(n-k) x (2n+1)] matrix;
% each row of Hin: first 2n bits give the Pauli, last entry is sign (1 or -1)

% Author: Narayanan Rengaswamy (July 1, 2021)


n = (size(Hin,2)-1)/2;
k = n - gfrank(Hin(:,1:2*n),2);

% Process the parity-check matrix to bring to the following form:
% H = [   0 , H_Z ; ... --> r_Z rows
%       H_A , H_B ]     --> r_X = (n-k)-r_Z rows
% All rows in this form must be linearly independent, and
% H_A must have full rank r_X

% First making sure that dependent rows are removed
Hin_full_rank = Hin;
if (size(Hin,1) > (n-k))
    redundant_rows = zeros(1,size(Hin,1) - (n-k));
    ind = 1;
    for i = 2:size(Hin,1)
        if (gfrank(Hin(1:i,1:2*n),2) == gfrank(Hin(1:(i-1),1:2*n),2))
            redundant_rows(ind) = i;
            ind = ind + 1;
        end
    end
    Hin_full_rank(redundant_rows,:) = [];
end

if (size(Hin_full_rank,1) ~= (n-k))
    fprintf('\nSome issue with rank of input H matrix!\n');
    exit;
end

% Process H s.t. the first r_Z rows are purely Z-type and the X part
% of the remaining r_X = (n-k)-r_Z rows has full rank r_X
r_X = gfrank(Hin_full_rank(:,1:n),2);
r_Z = (n-k) - r_X;
H = Hin_full_rank;
if (r_Z > 0)
    Zrows = zeros(1,r_Z);
    ind = 1;
    for i = 1:(n-k)
        if (all(H(i,1:n) == 0))
            Zrows(ind) = i;
            ind = ind + 1;
        elseif (i > 1) && (gfrank(H(1:i,1:n),2) == gfrank(H(1:(i-1),1:n),2))
            u = gflineq(H(1:(i-1),1:n)', H(i,1:n)', 2);
            H(i,:) = paulis_multiply(H(u==1,:));
            Zrows(ind) = i;
            ind = ind + 1;
        end
    end
    H = [ H(Zrows,:) ; H(setdiff(1:(n-k),Zrows),:) ];
end

% ------------------------------------------------------------------------------------------

% Simulate the creation of n Bell pairs
S_Bell = [ zeros(n), zeros(n),     eye(n),   eye(n),   ones(n,1) ; ...
            eye(n),   eye(n),   zeros(n), zeros(n),   ones(n,1) ];
      
% Simulate Alice's measurements      
replaced_rows = zeros(1,n-k);
Sout = S_Bell;
for i = 1:(n-k)
    [Sout, replaced_rows(i)] = stabilizer_formalism_msmt(Sout, [H(i,1:n), zeros(1,n),  H(i,(n+1):2*n), zeros(1,n),  H(i,end)]);
end

% Determine logical Z and X operators
Zbar = zeros(k,2*n);
inds = setdiff(1:n, replaced_rows);
logZ_ind = 1;
Hfull = H;
for i = 1:length(inds)
    Z_A = Sout(inds(i), (2*n+1):3*n);
    if (gfrank([Hfull(:,1:2*n); zeros(1,n), Z_A], 2) > gfrank(Hfull(:,1:2*n), 2))
        Zbar(logZ_ind, :) = [zeros(1,n), Z_A];
        Hfull(n-k+logZ_ind, :) = [zeros(1,n), Z_A,  1];
        logZ_ind = logZ_ind + 1;
    end
end

Xbar = zeros(k,2*n);
inds = setdiff((n+1):2*n, replaced_rows);
logX_ind = 1;
for i = 1:length(inds)
    X_A = [Sout(inds(i), 1:n) ,  Sout(inds(i), (2*n+1):3*n)];
    if (gfrank([Hfull(:,1:2*n); X_A], 2) > gfrank(Hfull(:,1:2*n), 2))
        Xbar(logX_ind, :) = X_A;
        Hfull(n+logX_ind, :) = [X_A,  1];
        logX_ind = logX_ind + 1;
    end
end

% Match the rows of Zbar accordingly with Xbar
symp_inn = mod(Zbar * fftshift(Xbar,2)', 2);
if (norm(symp_inn - eye(k),'fro') > 1e-10)
    symp_inn_inv = gf2matinv(symp_inn);
    Zbar = mod(symp_inn_inv * Zbar, 2);
end

end
