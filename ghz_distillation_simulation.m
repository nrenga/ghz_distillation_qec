% Script to perform one round of the GHZ distillation protocol using a
% given stabilizer code

% Prints the evolved GHZ stabilizers after each step of the protocol

% Author: Narayanan Rengaswamy (June 29, 2021)

clc
clear
% close all

%---------------------------------------------------------------------------------------------------------

% Choose if Alice (A) must do logical Clifford or Bob (B)
local_Clifford_by_A = 0;
local_Clifford_by_B = 1;

% Pauli Channel
epsilon = 0.1;
pvec = [ 1-epsilon, epsilon/3, epsilon/3, epsilon/3 ]; % [ p_I, p_X, p_Y, p_Z ]
pvec_cum = cumsum(pvec);

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

fprintf('\nGHZ DISTILLATION USING STABILIZER CODES\n');
fprintf('\n-------------------------------------------------\n');
fprintf('\nInput parity-check matrix for stabilizer code (last column is sign):\n');
print_matrix_with_spaces(Hin,n,2*n);
pause;

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
pause;

% [Xbar, Zbar] = find_logical_paulis_by_stab_msmt(H);
[Xbar, Zbar] = find_logical_paulis_by_ghz_msmt(H);

if ( any(any(mod(Zbar(:,1:2*n) * fftshift(Xbar(:,1:2*n),2)', 2) - eye(k))) )
    fprintf('\nSome issue with logical operators! Aborting ...\n\n');
    exit;
end

%---------------------------------------------------------------------------------------------------------

% Start with n copies of perfect GHZ states at Alice
S_GHZ = [ zeros(n), zeros(n), zeros(n),     eye(n),   eye(n), zeros(n),   ones(n,1) ; ...
          zeros(n), zeros(n), zeros(n),   zeros(n),   eye(n),   eye(n),   ones(n,1) ; ...
           eye(n),   eye(n),   eye(n),   zeros(n), zeros(n), zeros(n),   ones(n,1) ];

fprintf('\nStep 0: Stabilizers for n GHZ states (3n qubits) at Alice:\n');
print_matrix_with_spaces(S_GHZ, [1,2,4,5]*n, [3,6]*n, [n,2*n]);
pause;

% Simulate Alice's measurements
H_A = H;
replaced_rows_A = zeros(1,n-k);
syndrome_A = 1 - 2 * (rand(n-k,1)<0.5);
H_A(:,end) = syndrome_A;

fprintf('\nParity-check matrix for Alice''s code:\n');
print_matrix_with_spaces(H_A, [], [n,2*n], r_Z);

fprintf('\nLogical X operator(s):\n');
print_matrix_with_spaces(Xbar, [], [n,2*n]);

fprintf('\nLogical Z operator(s):\n');
print_matrix_with_spaces(Zbar, [] ,[n,2*n]);

pause;

S_A = S_GHZ;
for i = 1:(n-k)
    [S_A, replaced_rows_A(i)] = stabilizer_formalism_msmt(S_A, [H_A(i,1:n), zeros(1,2*n), ...
        H_A(i,(n+1):2*n), zeros(1,2*n),  H_A(i,end)]);
end

fprintf('\nStep 1: Alice measures her stabilizers on her qubits:\n');
print_matrix_with_spaces(S_A, [1,2,4,5]*n, [3,6]*n, [n,2*n]);
pause;

% Solve for binary symmetric matrix to map X components of stabilizers and
% logical X operators to their respective Z components; this ensures that
% Charlie also gets the same code, upto signs of stabilizers

% Standard symplectic matrix of the form
% T_R = [ I_n   R ;
%         0_n I_n ]
% R = R' (binary symmetric matrix)
%
% n = 3;
% R = [ 1 0 1 ;
%       0 1 0 ;
%       1 0 0 ]
% Ct_R = sqrt(Z)_1 * sqrt(Z)_2 * CZ_{13} --- Clifford circuit
%
% Ct_R * E(a,0) * Ct_R^\dagger = E(a, a*R) = E(a,b)
% a_1 = [1 1 0], b_1 = [1 1 0]; YYI
% a_2 = [0 1 1], b_2 = [0 1 1]; IYY
% a_3 = [1 0 0], b_3 = [1 0 0]; YII

A = [ H((r_Z+1):(n-k),1:n); Xbar(:,1:n) ];
B = [ H((r_Z+1):(n-k),(n+1):2*n); Xbar(:,(n+1):2*n) ];
R = find_binary_symmetric_matrix(A,B);

if (local_Clifford_by_A)
    R_A = R;  
else
    R_A = zeros(n);
end

T_R_A = [   eye(n) ,     R_A ; ...
        zeros(n) , eye(n) ];

S_A_Cliff = S_A;
S_A_Cliff(:,[(2*n+1):3*n, (5*n+1):6*n]) = mod(S_A(:,[(2*n+1):3*n, (5*n+1):6*n]) * T_R_A, 2);

if (local_Clifford_by_A)
    fprintf('\nStep 2: Alice applies diagonal Clifford to Charlie''s qubits:\n');
else
    fprintf('\nStep 2: Alice DOES NOT apply diagonal Clifford to Charlie''s qubits:\n');
end
print_matrix_with_spaces(S_A_Cliff, [1,2,4,5]*n, [3,6]*n, [n,2*n]);
pause;

%---------------------------------------------------------------------------------------------------------

% Construct BC stabilizers using the theorem on GHZ property
H_BC = [ H_A(1:r_Z,1:n) , zeros(r_Z,n) ,  H_A(1:r_Z,(n+1):2*n) , zeros(r_Z,n) ,  H_A(1:r_Z,end) ; ...
         H_A((r_Z+1):(n-k),1:n) , H_A((r_Z+1):(n-k),1:n)  ,  ...
         H_A((r_Z+1):(n-k),(n+1):2*n) , zeros(n-k-r_Z,n) ,  H_A((r_Z+1):(n-k),end) ; ...
          zeros(n) , zeros(n) ,  eye(n) , eye(n) ,  ones(n,1) ] ;
H_BC(:,[(n+1):2*n, (3*n+1):4*n]) = mod(H_BC(:,[(n+1):2*n, (3*n+1):4*n]) * T_R_A, 2);
signs_BC = [ H_A(:,end) .* (-1).^(sum(H_A(:,1:n) .* H_A(:,(n+1):2*n), 2)) ; ...
    ones(n,1) ] ;
H_BC(:,end) = signs_BC;

fprintf('\nParity-check matrix of joint code on BC qubits, induced by Alice:\n');
print_matrix_with_spaces(H_BC, [1,3]*n, [2,4]*n, [r_Z, n-k]);
pause;

% Sanity check: whether BC stabilizers are in S_A_Cliff
replaced_rows_BC = zeros(1,n-k);
S_BC = S_A_Cliff;
for i = 1:(n-k)
    [S_BC, replaced_rows_BC(i)] = stabilizer_formalism_msmt(S_BC, [zeros(1,n), H_BC(i,1:2*n), ...
        zeros(1,n), H_BC(i,(2*n+1):4*n),  H_BC(i,end)]);
end

if (any(any(S_BC - S_A_Cliff)))
    fprintf('\nProblem! Some BC stabilizers are not in the group!\n');
    exit;
end

%---------------------------------------------------------------------------------------------------------

% Simulate i.i.d. Pauli channel and Bob's 2n-qubit error correction
smpl = rand(1,2*n);
find_pauli = @(v) (find(v - pvec_cum' < 0, 1, 'first'));
err_BC = reshape(pauli(:,arrayfun(find_pauli,smpl))', 1, 4*n);

fprintf('\nStep 3: Channel introduces Pauli error on BC qubits:\n');
print_matrix_with_spaces(err_BC, [1,3]*n, [2,4]*n);

S_BC(:,end) = S_BC(:,end) .* (-1).^(mod(S_BC(:,[(n+1):3*n, (4*n+1):6*n]) * fftshift(err_BC,2)', 2));

fprintf('\nEffect of error on the stabilizers of all 3n qubits:\n');
print_matrix_with_spaces(S_BC, [1,2,4,5]*n, [3,6]*n, [n,2*n]);

syndrome_BC = mod(err_BC * fftshift(H_BC(:,1:4*n),2)', 2);

fprintf('\nSyndrome w.r.t. the joint BC code:\n');
print_matrix_with_spaces(syndrome_BC,[r_Z,n-k]);

pause;

err_est_BC = find_min_wt_error(H_BC, syndrome_BC);

fprintf('\nStep 4: Decoder estimates Pauli error on BC qubits:\n');
print_matrix_with_spaces(err_est_BC, [1,3]*n, [2,4]*n);

S_BC_corrected = S_BC;
S_BC_corrected(:,end) = S_BC_corrected(:,end) .* (-1).^(mod(S_BC_corrected(:,[(n+1):3*n, (4*n+1):6*n]) * fftshift(err_est_BC,2)', 2));

fprintf('\nEffect of error correction on the stabilizers of all 3n qubits:\n');
print_matrix_with_spaces(S_BC_corrected, [1,2,4,5]*n, [3,6]*n, [n,2*n]);
pause;

if (local_Clifford_by_B)
    R_B = R;
else
    R_B = zeros(n);
end

T_R_B = [   eye(n) ,     R_B ; ...
        zeros(n) , eye(n) ];

S_BC_Cliff = S_BC_corrected;
S_BC_Cliff(:,[(2*n+1):3*n, (5*n+1):6*n]) = mod(S_BC_corrected(:,[(2*n+1):3*n, (5*n+1):6*n]) * T_R_B, 2);

if (local_Clifford_by_B)
    fprintf('\nStep 5: Bob applies diagonal Clifford to Charlie''s qubits:\n');
else
    fprintf('\nStep 5: Bob DOES NOT apply diagonal Clifford to Charlie''s qubits:\n');
end
print_matrix_with_spaces(S_BC_Cliff, [1,2,4,5]*n, [3,6]*n, [n,2*n]);
pause;

% Simulate Bob's n-qubit measurements
H_B = H_A;
replaced_rows_B = zeros(1,n-k);
syndrome_B = 1 - 2 * (rand(n-k-r_Z,1)<0.5);
H_B((r_Z+1):(n-k),end) = H_A((r_Z+1):(n-k),end) .* syndrome_B;

S_B = S_BC_Cliff;
for i = 1:(n-k)
    [S_B, replaced_rows_B(i)] = stabilizer_formalism_msmt(S_B, [zeros(1,n), H_B(i,1:n), zeros(1,n), ...
        zeros(1,n), H_B(i,(n+1):2*n), zeros(1,n),  H_B(i,end)]);
end

fprintf('\nStep 6: Bob measures the code stabilizers on just his qubits:\n');
print_matrix_with_spaces(S_B, [1,2,4,5]*n, [3,6]*n, [n,2*n]);

fprintf('\nParity-check matrix of the code on B qubits:\n');
print_matrix_with_spaces(H_B, [], [n,2*n], r_Z);
pause;

%---------------------------------------------------------------------------------------------------------

% Simulate i.i.d. Pauli channel and Charlie's n-qubit error correction
smpl2 = rand(1,n);
err_C = reshape(pauli(:,arrayfun(find_pauli,smpl2))', 1, 2*n);

fprintf('\nStep 7: Channel introduces Pauli error on C qubits:\n');
print_matrix_with_spaces(err_C, [], [1,2]*n);

H_C = [ H_B(1:r_Z,:); H_B((r_Z+1):(n-k),1:n), zeros(n-k-r_Z,n),  H_B((r_Z+1):(n-k),end)];
H_C(:,1:2*n) = mod(H_C(:,1:2*n) * T_R_A * T_R_B, 2);
H_C((r_Z+1):(n-k),end) = H_BC((r_Z+1):(n-k),end) .* H_B((r_Z+1):(n-k),end);

XbarC = Xbar;
ZbarC = Zbar;
if (mod(local_Clifford_by_A + local_Clifford_by_B, 2) == 0)
    R_C = R;
    T_R_C = [   eye(n) ,     R_C ; ...
              zeros(n) , eye(n) ];
    XbarC(:,1:2*n) = mod(XbarC(:,1:2*n) * T_R_C, 2);
    XbarC(:,end) = XbarC(:,end) .* (-1).^(sum(XbarC(:,1:n) .* (XbarC(:,(n+1):2*n) .* (XbarC(:,1:n) * R_C)),2));
end

fprintf('\nParity-check matrix of the code on C qubits:\n');
print_matrix_with_spaces(H_C, [], [n,2*n], r_Z);

% Sanity check: whether C stabilizers are in S_B
replaced_rows_C = zeros(1,n-k);
S_C = S_B;
for i = 1:(n-k)
    [S_C, replaced_rows_C(i)] = stabilizer_formalism_msmt(S_C, [zeros(1,2*n), H_C(i,1:n), ...
        zeros(1,2*n), H_C(i,(n+1):2*n),  H_C(i,end)]);
end

if (any(any(S_C - S_B)))
    fprintf('\nProblem! Some C stabilizers are not in the group!\n');
    exit;
end

S_C(:,end) = S_C(:,end) .* (-1).^(mod(S_C(:,[(2*n+1):3*n, (5*n+1):6*n]) * fftshift(err_C,2)', 2));

fprintf('\nEffect of error on the stabilizers of all 3n qubits:\n');
print_matrix_with_spaces(S_C, [1,2,4,5]*n, [3,6]*n, [n,2*n]);

syndrome_C = mod(err_C * fftshift(H_C(:,1:2*n),2)', 2);

fprintf('\nSyndrome w.r.t. C code:\n');
print_matrix_with_spaces(syndrome_C,[r_Z,n-k]);

pause;

err_est_C = find_min_wt_error(H_C, syndrome_C);

fprintf('\nStep 8: Decoder estimates Pauli error on C qubits:\n');
print_matrix_with_spaces(err_est_C, [], [1,2]*n);

S_C_corrected = S_C;
S_C_corrected(:,end) = S_C_corrected(:,end) .* (-1).^(mod(S_C_corrected(:,[(2*n+1):3*n, (5*n+1):6*n]) * fftshift(err_est_C,2)', 2));

fprintf('\nEffect of error correction on the stabilizers of all 3n qubits:\n');
print_matrix_with_spaces(S_C_corrected, [1,2,4,5]*n, [3,6]*n, [n,2*n]);
pause;

%---------------------------------------------------------------------------------------------------------

% Check signs of logical GHZ stabilizers
GHZ_logical_ZZI = zeros(k, 6*n+1);
GHZ_logical_IZZ = zeros(k, 6*n+1);
GHZ_logical_XXX = zeros(k, 6*n+1);
GHZ_logical_ZZI_from_S_C = zeros(k, 6*n+1);
GHZ_logical_IZZ_from_S_C = zeros(k, 6*n+1);
GHZ_logical_XXX_from_S_C = zeros(k, 6*n+1);
v_ZZI = zeros(k, 3*n);
v_IZZ = zeros(k, 3*n);
v_XXX = zeros(k, 3*n);
for i = 1:k
    GHZ_logical_ZZI(i,:) = concatenate_pauli(concatenate_pauli(Zbar(i,:), Zbar(i,:)), [zeros(1,2*n), 1]);
    GHZ_logical_IZZ(i,:) = concatenate_pauli(concatenate_pauli([zeros(1,2*n), 1], Zbar(i,:)), ZbarC(i,:));
    GHZ_logical_XXX(i,:) = concatenate_pauli(concatenate_pauli(Xbar(i,:), Xbar(i,:)), XbarC(i,:));
    
    v_ZZI(i,:) = gflineq(S_C_corrected(:,1:6*n)', GHZ_logical_ZZI(i,1:6*n)')';
    GHZ_logical_ZZI_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_ZZI(i,:)==1,:));
    
    v_IZZ(i,:) = gflineq(S_C_corrected(:,1:6*n)', GHZ_logical_IZZ(i,1:6*n)')';
    GHZ_logical_IZZ_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_IZZ(i,:)==1,:));
    
    v_XXX(i,:) = gflineq(S_C_corrected(:,1:6*n)', GHZ_logical_XXX(i,1:6*n)')';
    GHZ_logical_XXX_from_S_C(i,:) = paulis_multiply(S_C_corrected(v_XXX(i,:)==1,:));
end

if ( any(any([GHZ_logical_ZZI(:,1:6*n); GHZ_logical_IZZ(:,1:6*n); GHZ_logical_XXX(:,1:6*n)] ...
        - [GHZ_logical_ZZI_from_S_C(:,1:6*n); GHZ_logical_IZZ_from_S_C(:,1:6*n); GHZ_logical_XXX_from_S_C(:,1:6*n)])) )
    fprintf('\nError! Some logical GHZ stabilizers were not produced correctly!\n');
    exit;
end

fprintf('-----------------------------------------------\n');
fprintf('\nREMINDER:\n');
fprintf('\nLogical X operator(s) for A and B:\n');
print_matrix_with_spaces(Xbar, [], [n,2*n]);
fprintf('\nLogical X operator(s) for C:\n');
print_matrix_with_spaces(XbarC, [], [n,2*n]);

fprintf('\nLogical Z operator(s) for A and B:\n');
print_matrix_with_spaces(Zbar, [], [n,2*n]);
fprintf('\nLogical Z operator(s) for C:\n');
print_matrix_with_spaces(ZbarC, [], [n,2*n]);
fprintf('-----------------------------------------------\n');

pause;

fprintf('\nThe logical GHZ stabilizer(s) of type ZZI exist in the group:\n');
print_matrix_with_spaces(GHZ_logical_ZZI_from_S_C, [1,2,4,5]*n, [3,6]*n);

fprintf('\nThe logical GHZ stabilizer(s) of type IZZ exist in the group:\n');
print_matrix_with_spaces(GHZ_logical_IZZ_from_S_C, [1,2,4,5]*n, [3,6]*n);

fprintf('\nThe logical GHZ stabilizer(s) of type XXX exist in the group:\n');
print_matrix_with_spaces(GHZ_logical_XXX_from_S_C, [1,2,4,5]*n, [3,6]*n);

pause;

%---------------------------------------------------------------------------------------------------------

S_decoded = S_C_corrected;

% Invert the Encoding Unitaries of Alice, Bob, and Charlie
sign_fix_A = fftshift(gflineq([Zbar(:,1:2*n); H_A(:,1:2*n); Xbar(:,1:2*n)], ...
    [Zbar(:,end); H_A(:,end); Xbar(:,end)]==-1), 2);
S_decoded(:,end) = S_decoded(:,end) .* (-1).^(mod(S_decoded(:,[1:n, (3*n+1):4*n]) * fftshift(sign_fix_A,2),2));
Dec_A = find_symp_mat([Zbar(:,1:2*n); H_A(:,1:2*n); Xbar(:,1:2*n)], ...
    [zeros(n), eye(n); eye(k), zeros(k,2*n-k)]);
S_decoded(:,[1:n, (3*n+1):4*n]) = mod(S_decoded(:,[1:n, (3*n+1):4*n]) * Dec_A, 2);

sign_fix_B = fftshift(gflineq([Zbar(:,1:2*n); H_B(:,1:2*n); Xbar(:,1:2*n)], ...
    [Zbar(:,end); H_B(:,end); Xbar(:,end)]==-1), 2);
S_decoded(:,end) = S_decoded(:,end) .* (-1).^(mod(S_decoded(:,[(n+1):2*n, (4*n+1):5*n]) * fftshift(sign_fix_B,2),2));
Dec_B = find_symp_mat([Zbar(:,1:2*n); H_B(:,1:2*n); Xbar(:,1:2*n)], ...
    [zeros(n), eye(n); eye(k), zeros(k,2*n-k)]);
S_decoded(:,[(n+1):2*n, (4*n+1):5*n]) = mod(S_decoded(:,[(n+1):2*n, (4*n+1):5*n]) * Dec_B, 2);

sign_fix_C = fftshift(gflineq([ZbarC(:,1:2*n); H_C(:,1:2*n); XbarC(:,1:2*n)], ...
    [ZbarC(:,end); H_C(:,end); XbarC(:,end)]==-1), 2);
S_decoded(:,end) = S_decoded(:,end) .* (-1).^(mod(S_decoded(:,[(2*n+1):3*n, (5*n+1):6*n]) * fftshift(sign_fix_C,2),2));
Dec_C = find_symp_mat([ZbarC(:,1:2*n); H_C(:,1:2*n); XbarC(:,1:2*n)], ...
    [zeros(n), eye(n); eye(k), zeros(k,2*n-k)]);
S_decoded(:,[(2*n+1):3*n, (5*n+1):6*n]) = mod(S_decoded(:,[(2*n+1):3*n, (5*n+1):6*n]) * Dec_C, 2);

fprintf('\nThe encoding unitary has been inverted separately on A''s, B''s, C''s qubits:\n');
print_matrix_with_spaces(S_decoded, [1,2,4,5]*n, [3,6]*n, [n,2*n]);

pause;

%---------------------------------------------------------------------------------------------------------

ZZI = [ zeros(k,3*n) ,  eye(k,n) , eye(k,n) , zeros(k,n) ,  ones(k,1) ];
IZZ = [ zeros(k,3*n) ,  zeros(k,n) , eye(k,n) , eye(k,n) ,  ones(k,1) ];
XXX = [ eye(k,n) , eye(k,n) , eye(k,n) ,  zeros(k,3*n) ,  ones(k,1) ];
ZZI_from_S_decoded = zeros(k, 6*n+1);
IZZ_from_S_decoded = zeros(k, 6*n+1);
XXX_from_S_decoded = zeros(k, 6*n+1);
w_ZZI = zeros(k, 3*n);
w_IZZ = zeros(k, 3*n);
w_XXX = zeros(k, 3*n);
for i = 1:k
    w_ZZI(i,:) = gflineq(S_decoded(:,1:6*n)', ZZI(i,1:6*n)')';
    ZZI_from_S_decoded(i,:) = paulis_multiply(S_decoded(w_ZZI(i,:)==1,:));
    
    w_IZZ(i,:) = gflineq(S_decoded(:,1:6*n)', IZZ(i,1:6*n)')';
    IZZ_from_S_decoded(i,:) = paulis_multiply(S_decoded(w_IZZ(i,:)==1,:));
    
    w_XXX(i,:) = gflineq(S_decoded(:,1:6*n)', XXX(i,1:6*n)')';
    XXX_from_S_decoded(i,:) = paulis_multiply(S_decoded(w_XXX(i,:)==1,:));
end

fprintf('\nThe GHZ stabilizer(s) of type ZZI exist in the group:\n');
print_matrix_with_spaces(ZZI_from_S_decoded, [1,2,4,5]*n, [3,6]*n);

fprintf('\nThe GHZ stabilizer(s) of type IZZ exist in the group:\n');
print_matrix_with_spaces(IZZ_from_S_decoded, [1,2,4,5]*n, [3,6]*n);

fprintf('\nThe GHZ stabilizer(s) of type XXX exist in the group:\n');
print_matrix_with_spaces(XXX_from_S_decoded, [1,2,4,5]*n, [3,6]*n);

pause;


% Check if the distillation was successful
% error_event = 0;
% if (any(GHZ_logical_ZZI(:,end) - GHZ_logical_ZZI_from_S_C(:,end)))
%     fprintf('\nDecoding error! Sign of a logical ZZI is incorrect.\n');
%     error_event = 1;
% end
% if (any(GHZ_logical_IZZ(:,end) - GHZ_logical_IZZ_from_S_C(:,end)))
%     fprintf('\nDecoding error! Sign of a logical IZZ is incorrect.\n');
%     error_event = 1;
% end
% if (any(GHZ_logical_XXX(:,end) - GHZ_logical_XXX_from_S_C(:,end)))
%     fprintf('\nDecoding error! Sign of a logical XXX is incorrect.\n');
%     error_event = 1;
% end

error_event = 0;
if (any(ZZI_from_S_decoded(:,end) == -1))
    fprintf('\nDecoding error! Sign of a ZZI is incorrect.\n');
    error_event = 1;
end
if (any(IZZ_from_S_decoded(:,end) == -1))
    fprintf('\nDecoding error! Sign of a IZZ is incorrect.\n');
    error_event = 1;
end
if (any(XXX_from_S_decoded(:,end) == -1))
    fprintf('\nDecoding error! Sign of a XXX is incorrect.\n');
    error_event = 1;
end

if (~error_event)
    fprintf('\nGHZ Distillation SUCCESSFUL!!\n\n');
else
    fprintf('\nGHZ Distillation FAILED!!\n\n');
end

