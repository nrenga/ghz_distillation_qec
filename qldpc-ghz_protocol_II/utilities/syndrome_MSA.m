function [decoder_output,Calculated_syndrome,decoding_success] = syndrome_MSA(H,syndrome,llr,max_iter)

[m,n] = size(H);

%% Initialization of channel probability.
vi = llr * ones(1,n);

%% Initialization  :  Assumption of a zero vector: channel LLR is sent from variable nodes to the check nodes in the zeroth iteration.
vij = llr * H;
cji = zeros(m,n);
lambda = 0.75;
    
for z = 1:max_iter
    %% Check node update :
    for i = 1:m
        row1s = find(H(i,:));
        row_wgt = length(row1s);
        row_vij = vij(i,row1s);
        for k = 1:row_wgt
            extrinsic_vij = row_vij;
            extrinsic_vij(k) = [];
            alpha = prod(sign(extrinsic_vij)) * (1-(2*syndrome(i)));
            cji(i,row1s(k)) = lambda * alpha * min(abs(extrinsic_vij));
        end
    end
    %% Variable node update:
    for j = 1:n
        col1s = find(cji(:,j));
        column_wgt = length(col1s);
        col_cji = cji(col1s,j);
        total_sum = sum(col_cji);
        Vi(j) =  vi(j) + total_sum;
        for k = 1:column_wgt
            extrinsic_cji = col_cji;
            extrinsic_cji(k) = [];
            vij(col1s(k),j) =  vi(j) + sum(extrinsic_cji);
        end
        decoder_output = (Vi < 0);   
    end
    
    %% Syndrome Calculation and stopping criteria
    Calculated_syndrome = (mod(H*(decoder_output)', 2));
    
    if all(Calculated_syndrome == syndrome)
        decoding_success = 1;
        break;
    else
        decoding_success = 0; %disp('Syndromes did not match!! Continue iteration');
    end
end % For iteration loop


end % For function.