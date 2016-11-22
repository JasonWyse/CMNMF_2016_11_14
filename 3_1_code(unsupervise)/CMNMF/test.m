function evaluation_test = test(matrix_cell_test,learned_matrix_cell,test_parameter_cell)
%TEST Summary of this function goes here
%   Detailed explanation goes here
%learned_matrix_cell = {L; W_out; H1_out; H2_out};
    W_out = learned_matrix_cell{2,1};    
    %matrix_cell_test =  {gene_idx;groundTruth_pathway_new;gene_idx_pathway;ppi_hsa_new;ppi_new_index;test_groundTruth_str};
    gene_idx = matrix_cell_test{1,1};
    groundTruth_pathway = matrix_cell_test{2,1};
    gene_idx_pathway = matrix_cell_test{3,1};
    ppi_new = matrix_cell_test{4,1};
    ppi_index = matrix_cell_test{5,1};
    test_groundTruth_str = matrix_cell_test{6,1};
    Z_threshold = 3;
    [W_Z_filtered,predicted_pathway_gene] = predicted_pathway(W_out,Z_threshold,gene_idx);
    %[Z_filter,pathway_ncbi]= predicted_pathway(W_out,T,mgi_id);
    F_measure_value = 1;
    predict_GG_matrix = W_Z_filtered*W_Z_filtered';
    if strcmp(test_groundTruth_str,'ppi') == 1
        [ F,jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_idx,ppi_new,...,
        ppi_index,F_measure_value);
    elseif strcmp(test_groundTruth_str,'pathway') == 1
        [ F,jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_idx,groundTruth_pathway,...,
        gene_idx_pathway,F_measure_value);
    end   
    M_sim = get_Ratio(predicted_pathway_gene);
    evaluation_test = [F,jaccard,RD,Precision,Recall,M_sim];

end

