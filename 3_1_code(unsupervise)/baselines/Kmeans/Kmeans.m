function [ Kmeans_result ] = Kmeans(species_acronym,file_cell,parameter_cell)
 
    g_p_network_file = file_cell{1,1};
    pathway_file = file_cell{2,1};
    ppi_file = file_cell{3,1};
    load(g_p_network_file,'g_p_network','gene_id');
    [m,n]=size(g_p_network);
    test_groundTruth_str = parameter_cell{1,1};
    K = parameter_cell{2,1};
    method_dir = parameter_cell{3,1};

    [Idx,~] = kmeans(g_p_network,K,'emptyaction','singleton');
    G_predict = zeros(m,K);
    for i = 1:m
        G_predict(i,Idx(i))=1; 
    end
    F_measure_value = 1;
    predict_GG_matrix = G_predict*G_predict';
    
    if strcmp(test_groundTruth_str,'pathway') ==1
        load(pathway_file, 'pathway_new');
        pathway_species_new =  pathway_new;
        [groundTruth_pathway, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_species_new);
        [F,Jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_id,groundTruth_pathway,...,
        gene_idx_pathway,F_measure_value);

        predicted_pathway_gene = transform(G_predict,gene_id);
        M_sim = get_Ratio(predicted_pathway_gene);
        Kmeans_result = [F,Jaccard,RD,Precision,Recall,M_sim];
        Kmeans_result_cell = {G_predict;Kmeans_result};
        
        result_file_name = ['Kmeans_filled_' species_acronym '_result_pathway_' datestr(now,30) '.mat' ];  
        save([method_dir result_file_name],'Kmeans_result_cell');
    end
    if strcmp(test_groundTruth_str,'ppi') ==1
        %load(ppi_file, 'ppi_new','ppi_new_index');
        load(ppi_file, 'ppi_difference','ppi_new_index');
        ppi_new = ppi_difference;
        [F,Jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_id,ppi_new,...,
        ppi_new_index,F_measure_value);
    
        predicted_pathway_gene = transform(G_predict,gene_id);
        M_sim = get_Ratio(predicted_pathway_gene);
        Kmeans_result = [F,Jaccard,RD,Precision,Recall,M_sim];
        Kmeans_result_cell = {G_predict;Kmeans_result};
        
        result_file_name = ['Kmeans_filled_' species_acronym '_result_ppi_difference' datestr(now,30) '.mat' ];  
        save([method_dir result_file_name],'Kmeans_result_cell');
    end
end

function[predicted_pathway_gene] = transform(learned_matrix,gene_idx)
    predicted_pathway_gene = zeros(size(learned_matrix,2),max(sum(learned_matrix)));
    for i = 1:size(predicted_pathway_gene,1)
       [row,~] = find(learned_matrix(:,i)==1); 
       predicted_pathway_gene(i,1:length(row)) = gene_idx(row);
    end
    predicted_pathway_gene = predicted_pathway_gene';
end