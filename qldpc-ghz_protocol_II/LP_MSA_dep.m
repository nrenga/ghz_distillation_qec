% function LP_MSA_dep(LPcode)

% Script to generate performance curves for a given
% Lifted Product code with the min-sum algorithm (MSA) decoder;
% These can be reinterpreted as the performance curves for Bell pair
% distillation depending on the noise setting as discussed in the paper
% https://arxiv.org/abs/2210.14143

% Author: Narayanan Rengaswamy (August 22, 2022)

clc
clear
close all

LPcode = 'LP_Matg8_dv3_dc5_L16_Dmin12_0'; % change according to the code

addpath('./utilities');
addpath('./codes');
addpath('./codes/LP_QC_QLDPC');

H_X = read_H(strcat(LPcode,'_hx'));
H_Z = read_H(strcat(LPcode,'_hz'));

H = [ zeros(size(H_Z)) , H_Z , ones(size(H_Z,1),1)  ; ...
      H_X , zeros(size(H_X)) , ones(size(H_X,1),1) ];

% [Xbar, Zbar] = find_logical_paulis_by_ghz_msmt(H);

L_X = read_L(strcat(LPcode,'_lx.alist'));
L_Z = read_L(strcat(LPcode,'_lz.alist'));

Xbar = [L_X, zeros(size(L_X))];
Zbar = [zeros(size(L_Z)), L_Z];

n = (size(H,2)-1)/2;
% k = n - gfrank(H(:,1:2*n),2);
k = size(L_X,1);

Omega_n = [ zeros(n) ,  eye(n) ; ...
             eye(n)  , zeros(n) ];

fprintf('\n[[%d, %d]] LP118 CODE + MSA DECODER\n',n,k);
fprintf('\n-------------------------------------------------\n');

if ( any(any(mod(Zbar(:,1:2*n) * Omega_n * Xbar(:,1:2*n)', 2) - eye(k))) )
    fprintf('\nSome issue with logical operators! Aborting ...\n\n');
    return;
end

% fprintf('\nLogical X operator(s):\n');
% print_matrix_with_spaces(Xbar, [], [n,2*n]);
% 
% fprintf('\nLogical Z operator(s):\n');
% print_matrix_with_spaces(Zbar, [] ,[n,2*n]);

% pause;

%---------------------------------------------------------------------------------------------------------

% epsilons = 0.1;
epsilons = [0.01; 0.03; 0.05; 0.07; 0.09; 0.1; 0.104; 0.108]
% epsilons = [1e-4; 5e-4; 1e-3; 4e-3; 7e-3; 0.01; 0.03; 0.05; 0.09 ];
% epsilons = [1e-5; 4e-5; 7e-5; 1e-4; 5e-4; 1e-3; 5e-3; 0.01; 0.06];
pvecs = [ 1-epsilons, epsilons/3, epsilons/3, epsilons/3 ]; % [ p_I, p_X, p_Y, p_Z ]
pvecs_cum = cumsum(pvecs,2);

% Number of simulations per data point, i.e., for each entry of 'epsilons'
blocks_min = 1e3
blocks_max = 1e6
no_of_logical_errors = 100

% Maximum number of iterations for MSA decoder
max_iter = 100

% Pauli:  I, X, Y, Z
pauli = [ 0, 1, 1, 0 ; ...
          0, 0, 1, 1 ];
      
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

fprintf('\n Starting simulation ... \n');

miscorrection = zeros(1,length(epsilons));
syndrome_unmatched = zeros(1,length(epsilons));

format longg

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
    
    while((iter < blocks_min) || ((miscorrection(eps_iter) + syndrome_unmatched(eps_iter)) < no_of_logical_errors))
        
        if (iter >= blocks_max)
            break;
        end
        
        iter = iter + 1;
    
        smpl = rand(1,n);
        err = reshape(pauli(:,arrayfun(find_pauli,smpl))', 1, 2*n);

        s_X = mod([zeros(size(H_Z)), H_Z] * err(1,[(n+1):2*n, 1:n])', 2);
        s_Z = mod([H_X, zeros(size(H_X))] * err(1,[(n+1):2*n, 1:n])', 2);
        
        [e_X_dec, s_X_dec, success_X_dec] = syndrome_MSA(H_Z, s_X, llr_X, max_iter);
        [e_Z_dec, s_Z_dec, success_Z_dec] = syndrome_MSA(H_X, s_Z, llr_Z, max_iter);
        
        if (success_X_dec && success_Z_dec)
            err_dec = [e_X_dec, e_Z_dec];
            err_eff = mod(err + err_dec, 2);
            if any(mod([Xbar(:,1:2*n) ; Zbar(:,1:2*n)] * err_eff(1,[(n+1):2*n, 1:n])',2))
                miscorrection(eps_iter) = miscorrection(eps_iter) + 1;
            end
        else
            syndrome_unmatched(eps_iter) = syndrome_unmatched(eps_iter) + 1;
        end
        
        if (mod(iter,250) == 0)
            %clc
            disp(epsilon);
            disp(iter);
            save(strcat('./LP_results/',LPcode,'_dep_MSA.mat'));
            miscorrection
            syndrome_unmatched
        end
        
    end
    miscorrection(eps_iter) = miscorrection(eps_iter)/iter;
    syndrome_unmatched(eps_iter) = syndrome_unmatched(eps_iter)/iter;
end

fprintf('\n Simulation complete! \n');

figure(1);
loglog(epsilons, miscorrection + syndrome_unmatched, '-ob');
xlabel('Depolarizing Probability', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Logical Error Rate', 'FontSize', 12, 'Interpreter', 'latex');
title('LP118 Code, Depolarizing Errors, MSA', 'FontSize', 12, 'Interpreter', 'latex');
grid on

% end    
    
   