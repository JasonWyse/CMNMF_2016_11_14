function [ output_args ] = NMF_ppi_difference()
%NMF_PPI_DIFFERENCE Summary of this function goes here
    species = 'human';
    if strcmp(species,'human') ==1
            species_acronym = 'hsa';         
    elseif strcmp(species,'mouse') == 1
        species_acronym = 'mmu';       
    end
    data_dir = ['../../../2_useful_data/' species '/'];
    tool_dir = '../../common_tool_functions/';
    path(path,data_dir); %invoke rmpath() to remove added path.
    path(path, tool_dir);
    load([data_dir 'g_p_network_2016_2.mat'],'g_p_network_first','g_p_network_second','M','gene_id');      
    load([data_dir [species_acronym '_pathway_data_2016_2.mat']], 'pathway_old');
    %load([data_dir [species_acronym '_pathway_data_2016_9.mat']], 'pathway_new');
    %load([data_dir [species_acronym '_ppi_2016_2.mat']], 'ppi_old','ppi_old_index'); 
    load([data_dir [species_acronym '_ppi_2016_0209_difference.mat']], 'ppi_difference','ppi_new_index');  

    dirOutput = dir('*.mat');
    fileNames={dirOutput.name}';
    name_prefix = ['NMF_' species_acronym '_result_ppi'];
    count =0;
    for i =1:length(fileNames)
       if(strncmp(fileNames(i),name_prefix,18) == 1)
           count = count+1;
           file_name_cell{count,1} = fileNames{i,1};
       end
    end
    load(file_name_cell{1,1},'NMF_result_cell');
    learned_matrix_cell = NMF_result_cell{2,1};
    best_parameter_array = NMF_result_cell{3,1};
    [dimension_K_set,~] = size(pathway_old);     
    max_ites = 80;%80 is a proper number
    top_n_set = [200;600;1000];
    test_groundTruth_str = 'ppi';
    matrix_cell_test =  {gene_id;ppi_difference;ppi_new_index;ppi_difference;ppi_new_index;test_groundTruth_str};
    test_parameter_cell = {best_parameter_array(1);best_parameter_array(2); dimension_K_set; max_ites ; top_n_set};
    evaluation_test = test(matrix_cell_test,learned_matrix_cell,test_parameter_cell);
    NMF_result_cell = {evaluation_test};
    %%%%%%%% save result %%%%%%%%
    result_file_name = ['NMF_' species_acronym '_result_ppi_difference_' datestr(now,30) '.mat' ];  
    save(result_file_name,'NMF_result_cell');
    %%%%%%%% save result %%%%%%%%
    rmpath(data_dir);
    rmpath(tool_dir);



end

