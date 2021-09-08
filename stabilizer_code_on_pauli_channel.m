% Script to simulate a stabilizer code on a standard Pauli channel

% Author: Narayanan Rengaswamy (July 27, 2021)

clc
clear
close all

%---------------------------------------------------------------------------------------------------------

% Pauli Channel 
epsilons = [1e-4; 5e-4; 1e-3; 5e-3; (0.01:0.02:0.1)'; 0.2];
pvecs = [ 1-epsilons, epsilons/3, epsilons/3, epsilons/3 ]; % [ p_I, p_X, p_Y, p_Z ]
pvecs_cum = cumsum(pvecs,2);

% Data for 5-qubit code
% epsilons = [1e-4; 5e-4; 1e-3; 5e-3; (0.01:0.02:0.1)'; 0.2];
% logical_error = [0, 0.0018, 0.0024, 0.0154, 0.0354, 0.0972, 0.1730,
%                     0.2272, 0.2950, 0.5664];

% Pauli:  I, X, Y, Z
pauli = [ 0, 1, 1, 0 ; ...
          0, 0, 1, 1 ];

%---------------------------------------------------------------------------------------------------------

% Stabilizer code parity-check matrix that has (2n+1) columns;
% last column is the sign for that row's stabilizer

% Hin = [ 0 0 0 ,  1 1 0 ,  1 ; ...
%         0 0 0 ,  0 1 1 ,  1 ];

% Hin = [ 1 1 0 ,  1 1 0 ,  1 ; ...
%         0 1 1 ,  0 1 1 ,  1 ];

% Hin = [ 1 1 0 ,  0 1 1 ,  1 ; ...
%         0 0 0 ,  1 1 0 ,  1 ; ...
%         1 1 0 ,  1 0 1 ,  1 ];

% Hin = [ 1 1 1 1 , 0 0 0 0 ,  1 ; ...
%         0 0 0 0 , 1 1 1 1 ,  1 ];

Hin = [ 1 0 0 1 0 ,  0 1 1 0 0 ,  1 ; ...
        0 1 0 0 1 ,  0 0 1 1 0 ,  1 ; ...
        1 0 1 0 0 ,  0 0 0 1 1 ,  1 ; ...
        0 1 0 1 0 ,  1 0 0 0 1 ,  1 ];

% Hin = [ 1 1 1 1 1 1 , 0 0 0 0 0 0 ,  1 ; ...
%         0 0 0 0 0 0 , 1 1 1 1 1 1 ,  1 ];

% Hin = [ 1 0 1 0 1 0 1 ,  0 0 0 0 0 0 0 ,  1 ; ...
%         0 1 1 0 0 1 1 ,  0 0 0 0 0 0 0 ,  1 ; ...
%         0 0 0 1 1 1 1 ,  0 0 0 0 0 0 0 ,  1 ; ...
%         0 0 0 0 0 0 0 ,  1 0 1 0 1 0 1 ,  1 ; ...
%         0 0 0 0 0 0 0 ,  0 1 1 0 0 1 1 ,  1 ; ...
%         0 0 0 0 0 0 0 ,  0 0 0 1 1 1 1 ,  1 ];

n = (size(Hin,2)-1)/2;
k = n - gfrank(Hin(:,1:2*n),2);

%---------------------------------------------------------------------------------------------------------

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

fprintf('\nParity-check matrix in standard form:\n');
print_matrix_with_spaces(H, [], [n,2*n], r_Z);
% pause;

% [Xbar, Zbar] = find_logical_paulis_by_stab_msmt(H);
[Xbar, Zbar] = find_logical_paulis_by_ghz_msmt(H);

if ( any(any(mod(Zbar(:,1:2*n) * fftshift(Xbar(:,1:2*n),2)', 2) - eye(k))) )
    fprintf('\nSome issue with logical operators! Aborting ...\n\n');
    exit;
end

fprintf('\nLogical X operator(s):\n');
print_matrix_with_spaces(Xbar, [], [n,2*n]);

fprintf('\nLogical Z operator(s):\n');
print_matrix_with_spaces(Zbar, [] ,[n,2*n]);

%---------------------------------------------------------------------------------------------------------

blocks = 5000;
logical_error = zeros(1,length(epsilons));
find_pauli = @(v) (find(v - pvec_cum' < 0, 1, 'first'));
        
for eps_iter = 1:length(epsilons)
    epsilon = epsilons(eps_iter)
    pvec = pvecs(eps_iter,:);
    pvec_cum = pvecs_cum(eps_iter,:);
    find_pauli = @(v) (find(v - pvec_cum' < 0, 1, 'first'));
    
    for iter = 1:blocks
        % Simulate i.i.d. Pauli channel and Charlie's n-qubit error correction
        smpl = rand(1,n);
        err = reshape(pauli(:,arrayfun(find_pauli,smpl))', 1, 2*n);
        
        syndrome = mod(err * fftshift(H(:,1:2*n),2)', 2);
        
        err_est = find_min_wt_error(H, syndrome);
        
        err_residue = mod(err + err_est, 2);
        
        if (any(mod(err_residue * fftshift([H(:,1:2*n); Xbar(:,1:2*n); Zbar(:,1:2*n)], 2)', 2)))
            logical_error(eps_iter) = logical_error(eps_iter) + 1;
        end
    end
    logical_error(eps_iter) = logical_error(eps_iter)/blocks
end

figure;
loglog(epsilons, logical_error, '-x'); 
xlabel('Depolarizing Probability ($\epsilon$)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Probability of Logical Error', 'FontSize', 12, 'Interpreter', 'latex');
title('Performance of the 5-Qubit Code', 'FontSize', 12, 'Interpreter', 'latex');
grid on

