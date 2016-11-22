function [ NMF_result_cell ] = NMF(parameter_cell,data_file_cell)
    g_p_network_file = data_file_cell{1,1};
    pathway_new_file = data_file_cell{2,1};
    ppi_new_file = data_file_cell{3,1};
    pathway_old_file = data_file_cell{4,1};
    ppi_old_file = data_file_cell{5,1};
    ppi_difference = data_file_cell{6,1};
    
    method_dir = parameter_cell{1,1};
    species = parameter_cell{2,1};    
    alpha_set = parameter_cell{3,1};
    beta_set  = parameter_cell{4,1};    
    max_ites = parameter_cell{5,1};
    top_n_set = parameter_cell{6,1};
    cv_criteria = parameter_cell{7,1};  
    validation_groundTruth_str = parameter_cell{8,1}; 
    test_groundTruth_str = validation_groundTruth_str;
    useful_data_dir = parameter_cell{9,1};
    baseline_name = parameter_cell{10,1};
    %%%%%%%%%%%%%%%%%%%%%%parse parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%load variable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    path(path,useful_data_dir);
    load(g_p_network_file ,'g_p_network_first','g_p_network_second','M','gene_id','first_level_id');      
    load(ppi_old_file, 'ppi_old','ppi_old_index'); 
    load(pathway_old_file, 'pathway_old');
    load(pathway_new_file, 'pathway_new');
    load(ppi_new_file,'ppi_new','ppi_new_index');  
    load(ppi_difference,'ppi_difference');
    %%%%%%%%%%%%%%%%%%%%%%load variable%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(species,'human')
        species_acronym = 'hsa';        
    elseif strcmp(species,'mouse')
        species_acronym = 'mmu';       
    end
    %%%%%%%%%%%%%%%%%%%%%% prepare data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V1 = g_p_network_first;
    V2 = g_p_network_second;    
    gene_idx = gene_id;
    [total_gene_num, total_p1_num] = size(V1);
    [~,total_p2_num] = size(V2);   
    [dimension_K_set,~] = size(pathway_old);   
    
    file_num_need = 10;
    CMNMF_data_dir = ['CMNMF/data/' species '/'];
    name_prefix = 'initial_matrix_fixed';
    get_files_parameter_cell = {file_num_need;CMNMF_data_dir;name_prefix;total_gene_num;...,
        dimension_K_set(1);total_p1_num;total_p2_num};
    file_name_cell_CMNMF = get_files(get_files_parameter_cell);
    CurrentMethod_data_dir = [method_dir 'data/' species '/'];
    Move_data_CMNMF_to_CurrentMethod(CMNMF_data_dir,CurrentMethod_data_dir);
    initial_matrixFileName_cell = file_name_cell_CMNMF; 
     %%%%%%%%%%%%%%%%%%%%%% prepare data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %%%%%%%%%%%%%%%%%%%%%%%%validation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [groundTruth_pathway_old, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_old);
    matrix_cell_validation  = {gene_idx;groundTruth_pathway_old;gene_idx_pathway;ppi_old;ppi_old_index;validation_groundTruth_str};
    input_parameter_cell = {alpha_set; beta_set; max_ites; top_n_set; cv_criteria; method_dir; species};    
    matrix_cell_train = {V1; V2; M; first_level_id};
    tic;
    [learned_matrix_cell,best_parameter_array,evaluation_result] = learn( input_parameter_cell, matrix_cell_train,...,
    matrix_cell_validation, initial_matrixFileName_cell );
    toc;
     %%%%%%%%%%%%%%%%%%%%%%%%validation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %%%%%%%%%%%%%%%test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [groundTruth_pathway_new, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_new);
    matrix_cell_test =  {gene_idx;groundTruth_pathway_new;gene_idx_pathway;ppi_new;ppi_new_index;test_groundTruth_str};
    test_parameter_cell = {best_parameter_array(1);best_parameter_array(2); dimension_K_set; max_ites ; top_n_set};
    evaluation_test = test(matrix_cell_test,learned_matrix_cell,test_parameter_cell);
    NMF_result_cell = {evaluation_test;learned_matrix_cell;best_parameter_array;evaluation_result;max_ites;alpha_set;beta_set};
    toc;
     %%%%%%%%%%%%%test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(test_groundTruth_str,'ppi') == 1
        result_file_name = [method_dir baseline_name '_' species_acronym '_result_ppi_' datestr(now,30) '.mat' ];  
        save(result_file_name, [baseline_name '_result_cell']);
    elseif strcmp(test_groundTruth_str,'pathway') == 1
        result_file_name = [method_dir baseline_name '_' species_acronym '_result_pathway_' datestr(now,30) '.mat' ];  
        save(result_file_name,[baseline_name '_result_cell']);
    end    
    rmpath(useful_data_dir); 

end
function [] = Move_data_CMNMF_to_CurrentMethod(CMNMF_data_dir,CurrentMethod_data_dir)
    if ~exist(CurrentMethod_data_dir,'dir')
        mkdir(CurrentMethod_data_dir);        
    end
    copyfile([CMNMF_data_dir '*'],CurrentMethod_data_dir);
    
end
function [file_name_cell] = get_files(get_files_parameter_cell)    
    file_num_need = get_files_parameter_cell{1,1};
    dir_path = get_files_parameter_cell{2,1};
    name_prefix = get_files_parameter_cell{3,1};
    total_gene_num = get_files_parameter_cell{4,1};
    dimension_K_set(1) = get_files_parameter_cell{5,1};
    total_p1_num = get_files_parameter_cell{6,1};
    total_p2_num = get_files_parameter_cell{7,1};
    
    fileFolder=fullfile(dir_path);
    dirOutput=dir(fullfile(fileFolder,'*.mat'));
    fileNames={dirOutput.name}';
    count =0;
    for i =1:length(fileNames)
       if(strncmp(fileNames(i),name_prefix,5) == 1)
           count = count+1;
           file_name_cell{count,1} = fileNames{i,1};
       end
    end
    if count > file_num_need
       file_name_cell_temp = file_name_cell(end-file_num_need+1:end,:); 
       file_name_cell = file_name_cell_temp;
    elseif count < file_num_need
        count_temp = count;
        for i = 1:file_num_need - count
            file_name_new = ['data/initial_matrix_fixed_' datestr(now,30) '.mat' ];            
            W = rand(total_gene_num, dimension_K_set(1));
            H1 = rand(dimension_K_set(1), total_p1_num);
            H2 = rand(dimension_K_set(1), total_p2_num);    
            save(file_name_new,'W','H1','H2');
            count_temp = count_temp + 1;
            file_name_cell{count_temp,1} = file_name_new;
        end
    end
end