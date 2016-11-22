function [ CMNMF_result_cell ] = CMNMF_main()
%CMNMF_MAIN Summary of this function goes here
%   Detailed explanation goes here   
    method_dir = 'CMNMF/';
    species = 'mouse';    
    validation_groundTruth_str = 'ppi';%two choices: 'pathway' or 'ppi'
    alpha_set = [0.001,0.01,0.1,1,10,100,1000];%0.001,0.01,0.1,1,10,100,1000
    beta_set  = [0.001,0.01,0.1,1,10,100,1000];%0.001,0.01,0.1,1,10,100,1000
    %dimension_K_set = [20];
    max_ites = 80;%80 is a proper number
    top_n_set = [200;600;1000];
    cv_criteria = 'F1';    
    parameter_cell = {method_dir;species;alpha_set;
                        beta_set;max_ites;top_n_set;
                        cv_criteria;validation_groundTruth_str};
    path(path,method_dir);
    path(path,'common_tool_functions');    
    CMNMF_result_cell = CMNMF(parameter_cell);
    validation_groundTruth_str = 'pathway';
    parameter_cell{8,1} = validation_groundTruth_str;
    CMNMF_result_cell = CMNMF(parameter_cell);
    rmpath(method_dir);
    rmpath('common_tool_functions');
end

