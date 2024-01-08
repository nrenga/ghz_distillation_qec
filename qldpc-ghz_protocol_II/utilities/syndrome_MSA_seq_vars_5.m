function [decoder_output,Calculated_syndrome,decoding_success] = syndrome_MSA_seq_vars_5(H,vn_neighbors,cn_neighbors,syndrome,llr,max_iter)

[m,n] = size(H);
syndromeBPSK = 1-2*syndrome;

%% Initialization of channel probability.
vj = llr * ones(1,n);
Vj = zeros(1,n);

%% Initialization  :  Assumption of a zero vector: channel LLR is sent from variable nodes to the check nodes in the zeroth iteration.
vji = llr * H;
cij = zeros(m,n);
lambda = 0.8;

%% Check node update:
for i = 1:m

    row1s = cn_neighbors{i};
    row_vji = vji(i,row1s);
%     extrinsic_vjis = toeplitz([row_vji(end) row_vji(1:(end-1))], fliplr(row_vji));
%     extrinsic_vjis(:,end) = [];
    extrinsic_vjis = (ones(length(row1s),1) * row_vji) .* (1-eye(length(row1s))); 
    extrinsic_vjis = extrinsic_vjis + 1e4*eye(length(row1s));    
    alphas = syndromeBPSK(i) * prod(sign(extrinsic_vjis),2)';
    betas = min(abs(extrinsic_vjis),[],2)';
    cij(i,row1s) = lambda * (alphas .* betas);

end

    
for z = 1:max_iter
    %% Variable node update:
    for j = 1:n

        col1s = vn_neighbors{j};
        
        for i = 1:length(col1s)

            col_cij_tot = sum(cij(col1s,j));
            vji(col1s(i),j) = vj(j) + col_cij_tot - cij(col1s(i),j);

            %% Check node update:
            row1s = cn_neighbors{col1s(i)};
            row_vji = vji(col1s(i),row1s);
            %extrinsic_vjis = toeplitz([row_vji(end) row_vji(1:(end-1))], fliplr(row_vji));
            %extrinsic_vjis = [circshift(row_vji,-1); zeros(length(row1s)-1,length(row1s))];
            %L = circshift(eye(length(row1s)),1);
            %for k = 1:length(row1s)
            %    extrinsic_vjis(k,:) = row_vji * (L^k);
                %extrinsic_vjis(k,:) = circshift(extrinsic_vjis(k-1,:),-1);
            %end
            %extrinsic_vjis(:,end) = [];
            
            extrinsic_vjis = (ones(length(row1s),1) * row_vji) .* (1-eye(length(row1s)));
            extrinsic_vjis = extrinsic_vjis + 1e4*eye(length(row1s));
            alphas = syndromeBPSK(col1s(i)) * prod(sign(extrinsic_vjis),2)';
            betas = min(abs(extrinsic_vjis),[],2)';
            cij(col1s(i),row1s) = lambda * (alphas .* betas);
            
        end
        
        %% Variable node decisions:
        total_sum = sum(cij(col1s,j));
        Vj(j) =  vj(j) + total_sum;
        decoder_output = (Vj < 0);
        
        %% Syndrome Calculation and stopping criteria
        Calculated_syndrome = mod(H*(decoder_output)', 2);
        
        if all(Calculated_syndrome == syndrome)
            decoding_success = 1;
            return;
        else
            decoding_success = 0; %disp('Syndromes did not match!! Continue iteration');
        end
        
    end
end % For iteration loop


end % For function.