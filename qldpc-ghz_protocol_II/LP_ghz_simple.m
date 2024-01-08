% function LP_ghz_simple(LPcode)

% Script to generate performance curves for GHZ distillation using a given
% Lifted Product code and the simpler Protocol II in the paper
% https://arxiv.org/abs/2210.14143

% Author: Narayanan Rengaswamy (August 22, 2022)

clc
clear
close all

LPcode = 'LP_Matg8_dv3_dc5_L16_Dmin12_0'; % change according to the code

addpath('./utilities');
addpath('./codes');
addpath('./codes/LP_QC_QLDPC');

%---------------------------------------------------------------------------------------------------------

% Pauli Channel
% epsilons = [1e-4; 5e-4; 1e-3; 5e-3; (0.01:0.02:0.1)'; 0.2];
% epsilons = [(1:9)'*1e-4; 3e-3; 7e-3; 0.01 ];
epsilons = [0.05; 0.07; 0.09; 0.1; 0.104; 0.108; 0.11];

% epsilons = [0.05; 0.07; 0.09; 0.1; 0.104; 0.108; 0.11];
% epsilons = [ 1e-4; 5e-4; 1e-3; 4e-3; 7e-3; 0.01; 0.03; 0.05; 0.09 ];

epsilons
pvecs = [ 1-epsilons, epsilons/3, epsilons/3, epsilons/3 ]; % [ p_I, p_X, p_Y, p_Z ]
pvecs_cum = cumsum(pvecs,2);

% Number of simulations per data point, i.e., for each entry of 'epsilons'
blocks_min = 1e3
blocks_max = 1e5
no_of_logical_errors = 10000

% Maximum number of iterations for MSA decoder
max_iter = 100

% Pauli:  I, X, Y, Z
pauli = [ 0, 1, 1, 0 ; ...
          0, 0, 1, 1 ];

%---------------------------------------------------------------------------------------------------------

% Stabilizer code parity-check matrix that has (2n+1) columns;
% last column is the sign for that row's stabilizer

H_X = read_H(strcat(LPcode,'_hx'));
H_Z = read_H(strcat(LPcode,'_hz'));

H = [ zeros(size(H_Z)) , H_Z , ones(size(H_Z,1),1)  ; ...
      H_X , zeros(size(H_X)) , ones(size(H_X,1),1) ];

% [Xbar, Zbar] = find_logical_paulis_by_ghz_msmt(H);

L_X = read_L(strcat(LPcode,'_lx.alist'));
L_Z = read_L(strcat(LPcode,'_lz.alist'));

n = (size(H,2)-1)/2;
% k = n - gfrank(H(:,1:2*n),2);
k = size(L_X,1);
r = size(H,1);
r_Z = size(H_Z,1);

Xbar = [L_X, zeros(size(L_X)), ones(k,1)];
Zbar = [zeros(size(L_Z)), L_Z, ones(k,1)];

Omega_n = [ zeros(n) ,  eye(n) ; ...
             eye(n)  , zeros(n) ];

fprintf('\nSIMPLER GHZ DISTILLATION: [[%d, %d]] LP118 CODE + MSA DECODER\n',n,k);
fprintf('\n-------------------------------------------------\n');

if ( any(any(mod(Zbar(:,1:2*n) * Omega_n * Xbar(:,1:2*n)', 2) - eye(k))) )
    fprintf('\nSome issue with logical operators! Aborting ...\n\n');
    return;
end

%---------------------------------------------------------------------------------------------------------

H_X_vn_neighbors = cell(n,1);
H_Z_vn_neighbors = cell(n,1);
H_X_cn_neighbors = cell(size(H_X,1),1);
H_Z_cn_neighbors = cell(size(H_Z,1),1);

for j = 1:n
    H_X_vn_neighbors{j} = find(H_X(:,j));
    H_Z_vn_neighbors{j} = find(H_Z(:,j));
end

for i = 1:max(size(H_X,1),size(H_Z,1))
    if (i <= size(H_X,1))
        H_X_cn_neighbors{i} = find(H_X(i,:))';
    end
    
    if (i <= size(H_Z,1))
        H_Z_cn_neighbors{i} = find(H_Z(i,:))';
    end
end

logical_error_ghz = zeros(1,length(epsilons));

GHZ_logical_ZZI = zeros(k, 6*n+1);
GHZ_logical_IZZ = zeros(k, 6*n+1);
GHZ_logical_XXX = zeros(k, 6*n+1);

for i = 1:k
    GHZ_logical_ZZI(i,:) = concatenate_pauli(concatenate_pauli(Zbar(i,:), Zbar(i,:)), [zeros(1,2*n), 1]);
    GHZ_logical_IZZ(i,:) = concatenate_pauli(concatenate_pauli([zeros(1,2*n), 1], Zbar(i,:)), Zbar(i,:));
    GHZ_logical_XXX(i,:) = concatenate_pauli(concatenate_pauli(Xbar(i,:), Xbar(i,:)), Xbar(i,:));
end    

v_ZZI = zeros(k, 3*n);
v_IZZ = zeros(k, 3*n);
v_XXX = zeros(k, 3*n);

%---------------------------------------------------------------------------------------------------------

% Simulate one round of distillation to obtain all relevant information

% Start with n copies of perfect GHZ states at Alice
S_GHZ = [ zeros(n), zeros(n), zeros(n),     eye(n),   eye(n), zeros(n),   ones(n,1) ; ...
          zeros(n), zeros(n), zeros(n),   zeros(n),   eye(n),   eye(n),   ones(n,1) ; ...
           eye(n),   eye(n),   eye(n),   zeros(n), zeros(n), zeros(n),   ones(n,1) ];

% Simulate Alice's measurements on her qubits
H_A = H;
replaced_rows_A = zeros(1,r);
syndrome_A = 1 - 2 * (rand(r,1)<0.5);
H_A(:,end) = syndrome_A;

S_A = S_GHZ;
for i = 1:r
    [S_A, replaced_rows_A(i)] = stabilizer_formalism_msmt(S_A, [H_A(i,1:n), zeros(1,2*n), ...
        H_A(i,(n+1):2*n), zeros(1,2*n),  H_A(i,end)]);
end

%---------------------------------------------------------------------------------------------------------

% Construct BC stabilizers using the theorem on GHZ property
H_BC = [ H_A(1:r_Z,1:n) , zeros(r_Z,n) ,  H_A(1:r_Z,(n+1):2*n) , zeros(r_Z,n) ,  H_A(1:r_Z,end) ; ...
         H_A((r_Z+1):r,1:n) , H_A((r_Z+1):r,1:n)  ,  ...
         zeros(r-r_Z,n) , zeros(r-r_Z,n) ,  H_A((r_Z+1):r,end) ; ...
         zeros(n) , zeros(n) ,  eye(n) , eye(n) ,  ones(n,1) ] ;

%---------------------------------------------------------------------------------------------------------

% Simulate Alice's measurements on Bob's qubits
H_B = H_A;
replaced_rows_B = zeros(1,r);
syndrome_AB = 1 - 2 * (rand(r-r_Z,1)<0.5);
H_B((r_Z+1):r,end) = H_A((r_Z+1):r,end) .* syndrome_AB;

S_B = S_A;
for i = 1:r
    [S_B, replaced_rows_B(i)] = stabilizer_formalism_msmt(S_B, [zeros(1,n), H_B(i,1:n), zeros(1,n), ...
        zeros(1,n), H_B(i,(n+1):2*n), zeros(1,n),  H_B(i,end)]);
end

%---------------------------------------------------------------------------------------------------------

% Simulate i.i.d. Pauli channel and Bob's n-qubit error correction
err_B = zeros(1,2*n);

S_B(:,end) = S_B(:,end) .* (-1).^(mod(S_B(:,[(4*n+1):5*n, (n+1):2*n]) * err_B', 2));

syndrome_B = mod(err_B * H_B(:,[(n+1):2*n, 1:n])', 2);

errX_est_B = syndrome_MSA_seq_vars_5(H_B(1:r_Z,(n+1):2*n), syndrome_B(1:r_Z), 0, max_iter);
errZ_est_B = syndrome_MSA_seq_vars_5(H_B((r_Z+1):r,1:n), syndrome_B((r_Z+1):r), 0, max_iter);
err_est_B = [errX_est_B, errZ_est_B];

S_B_corrected = S_B;
S_B_corrected(:,end) = S_B_corrected(:,end) .* (-1).^(mod(S_B_corrected(:,[(4*n+1):5*n, (n+1):2*n]) * err_est_B', 2));

%---------------------------------------------------------------------------------------------------------

% Simulate i.i.d. Pauli channel and Charlie's n-qubit error correction
err_C = zeros(1,2*n);

H_C = H_B;
H_C((r_Z+1):r,end) = H_BC((r_Z+1):r,end) .* H_B((r_Z+1):r,end);

S_C = S_B_corrected;
S_C(:,end) = S_C(:,end) .* (-1).^(mod(S_C(:,[(5*n+1):6*n, (2*n+1):3*n]) * err_C', 2));

syndrome_C = mod(err_C * H_C(:,[(n+1):2*n, 1:n])', 2);

errX_est_C = syndrome_MSA_seq_vars_5(H_C(1:r_Z,(n+1):2*n), syndrome_C(1:r_Z), 0, max_iter);
errZ_est_C = syndrome_MSA_seq_vars_5(H_C((r_Z+1):r,1:n), syndrome_C((r_Z+1):r), 0, max_iter);
err_est_C = [errX_est_C, errZ_est_C];

S_C_corrected = S_C;
S_C_corrected(:,end) = S_C_corrected(:,end) .* (-1).^(mod(S_C_corrected(:,[(5*n+1):6*n, (2*n+1):3*n]) * err_est_C', 2));

%---------------------------------------------------------------------------------------------------------

% Check signs of logical GHZ stabilizers
GHZ_logical_ZZI_from_S_C = zeros(k, 6*n+1);
GHZ_logical_IZZ_from_S_C = zeros(k, 6*n+1);
GHZ_logical_XXX_from_S_C = zeros(k, 6*n+1);
for i = 1:k
    v_ZZI(i,:) = [gflineq(S_C_corrected(1:n,1:6*n)', GHZ_logical_ZZI(i,1:6*n)')' , zeros(1,2*n)];
    v_IZZ(i,:) = [zeros(1,n) , gflineq(S_C_corrected((n+1):2*n,1:6*n)', GHZ_logical_IZZ(i,1:6*n)')' , zeros(1,n)];
    v_XXX(i,:) = [zeros(1,2*n) , gflineq(S_C_corrected((2*n+1):3*n,1:6*n)', GHZ_logical_XXX(i,1:6*n)')'];

    GHZ_logical_ZZI_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_ZZI(i,:)==1,:));
    GHZ_logical_IZZ_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_IZZ(i,:)==1,:));
    GHZ_logical_XXX_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_XXX(i,:)==1,:));
end

if ( any(any([GHZ_logical_ZZI(:,1:6*n); GHZ_logical_IZZ(:,1:6*n); GHZ_logical_XXX(:,1:6*n)] ...
        - [GHZ_logical_ZZI_from_S_C(:,1:6*n); GHZ_logical_IZZ_from_S_C(:,1:6*n); GHZ_logical_XXX_from_S_C(:,1:6*n)])) )
    fprintf('\nError! Some logical GHZ stabilizers were not produced correctly!\n');
    return;
end

S_B_in = S_B;

%---------------------------------------------------------------------------------------------------------


fprintf('\n Starting simulation ... \n');

format longg

% pause;

for eps_iter = 1:length(epsilons)
    
    epsilon = epsilons(eps_iter);
    pvec = pvecs(eps_iter,:);
    
    p_Z = pvec(3) + pvec(4);
    llr_Z = log((1-p_Z)/p_Z);
    p_X = pvec(2) + pvec(3);
    llr_X = log((1-p_X)/p_X);
    
    pvec_cum = pvecs_cum(eps_iter,:);
    find_pauli = @(v) (find(v - pvec_cum' < 0, 1, 'first'));
    
    iter = 0;
    
    while((iter < blocks_min) || (logical_error_ghz(eps_iter) < no_of_logical_errors))
        
        if (iter >= blocks_max)
            break;
        end
        
        iter = iter + 1;

        %---------------------------------------------------------------------------------------------------------
        
        % Start with n copies of perfect GHZ states at Alice
%         S_GHZ = [ zeros(n), zeros(n), zeros(n),     eye(n),   eye(n), zeros(n),   ones(n,1) ; ...
%                   zeros(n), zeros(n), zeros(n),   zeros(n),   eye(n),   eye(n),   ones(n,1) ; ...
%                     eye(n),   eye(n),   eye(n),   zeros(n), zeros(n), zeros(n),   ones(n,1) ];
%         
        % Simulate Alice's measurements on her qubits
        H_A = H;
        syndrome_A = 1 - 2 * (rand(r,1)<0.5);
        H_A(:,end) = syndrome_A;
        
%         for i = 1:r
%             S_A(replaced_rows_A(i),:) = [H_A(i,1:n), zeros(1,2*n), H_A(i,(n+1):2*n), zeros(1,2*n),  H_A(i,end)];
%         end
%         
        %---------------------------------------------------------------------------------------------------------
        
        % Construct BC stabilizers using the theorem on GHZ property
        H_BC = [ H_A(1:r_Z,1:n) , zeros(r_Z,n) ,  H_A(1:r_Z,(n+1):2*n) , zeros(r_Z,n) ,  H_A(1:r_Z,end) ; ...
                 H_A((r_Z+1):r,1:n) , H_A((r_Z+1):r,1:n)  ,  ...
                 zeros(r-r_Z,n) , zeros(r-r_Z,n) ,  H_A((r_Z+1):r,end) ; ...
                 zeros(n) , zeros(n) ,  eye(n) , eye(n) ,  ones(n,1) ] ;

        %---------------------------------------------------------------------------------------------------------

        % Simulate Alice's measurements on Bob's qubits
        H_B = H_A;
        syndrome_AB = 1 - 2 * (rand(r-r_Z,1)<0.5);
        H_B((r_Z+1):r,end) = H_A((r_Z+1):r,end) .* syndrome_AB;

        S_B = S_B_in;
        for i = 1:r
            if (replaced_rows_A(i) > 0)
                S_B(replaced_rows_A(i),:) = [H_A(i,1:n), zeros(1,2*n), H_A(i,(n+1):2*n), zeros(1,2*n),  H_A(i,end)];
            end
            if (replaced_rows_B(i) > 0)
                S_B(replaced_rows_B(i),:) = [zeros(1,n), H_B(i,1:n), zeros(1,n), zeros(1,n), H_B(i,(n+1):2*n), zeros(1,n),  H_B(i,end)];
            end
        end
        
        %---------------------------------------------------------------------------------------------------------
        
        % Simulate i.i.d. Pauli channel and Bob's n-qubit error correction
        smpl = rand(1,n);
        err_B = reshape(pauli(:,arrayfun(find_pauli,smpl))', 1, 2*n);
        
        S_B(:,end) = S_B(:,end) .* (-1).^(mod(S_B(:,[(4*n+1):5*n, (n+1):2*n]) * err_B', 2));
        
        syndrome_B = mod(err_B * H_B(:,[(n+1):2*n, 1:n])', 2);
        
        errX_est_B = syndrome_MSA_seq_vars_5(H_B(1:r_Z,(n+1):2*n), syndrome_B(1:r_Z), llr_X, max_iter);
        errZ_est_B = syndrome_MSA_seq_vars_5(H_B((r_Z+1):r,1:n), syndrome_B((r_Z+1):r), llr_Z, max_iter);
        err_est_B = [errX_est_B, errZ_est_B];
        
        S_B_corrected = S_B;
        S_B_corrected(:,end) = S_B_corrected(:,end) .* (-1).^(mod(S_B_corrected(:,[(4*n+1):5*n, (n+1):2*n]) * err_est_B', 2));
        
        %---------------------------------------------------------------------------------------------------------

        % Simulate i.i.d. Pauli channel and Charlie's n-qubit error correction
        smpl2 = rand(1,n);
        err_C = reshape(pauli(:,arrayfun(find_pauli,smpl2))', 1, 2*n);
        
        H_C = H_B;
        H_C((r_Z+1):r,end) = H_BC((r_Z+1):r,end) .* H_B((r_Z+1):r,end);
        
        S_C = S_B_corrected;
        S_C(:,end) = S_C(:,end) .* (-1).^(mod(S_C(:,[(5*n+1):6*n, (2*n+1):3*n]) * err_C', 2));
        
        syndrome_C = mod(err_C * H_C(:,[(n+1):2*n, 1:n])', 2);
        
        errX_est_C = syndrome_MSA_seq_vars_5(H_C(1:r_Z,(n+1):2*n), syndrome_C(1:r_Z), llr_X, max_iter);
        errZ_est_C = syndrome_MSA_seq_vars_5(H_C((r_Z+1):r,1:n), syndrome_C((r_Z+1):r), llr_Z, max_iter);
        err_est_C = [errX_est_C, errZ_est_C];
        
        S_C_corrected = S_C;
        S_C_corrected(:,end) = S_C_corrected(:,end) .* (-1).^(mod(S_C_corrected(:,[(5*n+1):6*n, (2*n+1):3*n]) * err_est_C', 2));
        
        %---------------------------------------------------------------------------------------------------------
        
%         % Check signs of logical GHZ stabilizers
        GHZ_logical_ZZI_from_S_C = zeros(k, 6*n+1);
        GHZ_logical_IZZ_from_S_C = zeros(k, 6*n+1);
        GHZ_logical_XXX_from_S_C = zeros(k, 6*n+1);
        for i = 1:k
            GHZ_logical_ZZI_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_ZZI(i,:)==1,:));
            GHZ_logical_IZZ_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_IZZ(i,:)==1,:));
            GHZ_logical_XXX_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_XXX(i,:)==1,:));
        end
        
        if ( any(any([GHZ_logical_ZZI(:,1:6*n); GHZ_logical_IZZ(:,1:6*n); GHZ_logical_XXX(:,1:6*n)] ...
                - [GHZ_logical_ZZI_from_S_C(:,1:6*n); GHZ_logical_IZZ_from_S_C(:,1:6*n); GHZ_logical_XXX_from_S_C(:,1:6*n)])) )
            fprintf('\nError! Some logical GHZ stabilizers were not produced correctly!\n');
            return;
        end

        %---------------------------------------------------------------------------------------------------------
        
        % Check if the distillation was successful
        error_event = 0;
        if (any(GHZ_logical_ZZI(:,end) - GHZ_logical_ZZI_from_S_C(:,end)))
%             fprintf('\nDecoding error! Sign of a logical ZZI is incorrect.\n');
            error_event = 1;
        end
        if (any(GHZ_logical_IZZ(:,end) - GHZ_logical_IZZ_from_S_C(:,end)))
%             fprintf('\nDecoding error! Sign of a logical IZZ is incorrect.\n');
            error_event = 1;
        end
        if (any(GHZ_logical_XXX(:,end) - GHZ_logical_XXX_from_S_C(:,end)))
%             fprintf('\nDecoding error! Sign of a logical XXX is incorrect.\n');
            error_event = 1;
        end
        
        logical_error_ghz(eps_iter) = logical_error_ghz(eps_iter) + error_event;
        
        if (mod(iter,25) == 0)
            clc
            disp(epsilon);
            disp(iter);
            save(strcat('./LP_results/GHZsimple_',LPcode,'_dep_MSA.mat'));
            logical_error_ghz
        end
        
    end

    logical_error_ghz(eps_iter) = logical_error_ghz(eps_iter)/iter;
    
end

save(strcat('./LP_results/GHZsimple_',LPcode,'_dep_MSA.mat'));
fprintf('\n Simulation complete! \n');

figure;
loglog(epsilons, logical_error_ghz, '-ob');
xlabel('Depolarizing Probability', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Probability of Failure', 'FontSize', 12, 'Interpreter', 'latex');
title('GHZ Distillation Protocol II with the LP118 Code, Depolarizing Errors, MSA', 'FontSize', 12, 'Interpreter', 'latex');
% xlim([1e-4,1e-1]);
% ylim([1e-4,1]);
grid on

% end
