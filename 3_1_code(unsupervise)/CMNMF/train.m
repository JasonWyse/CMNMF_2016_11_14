function evaluation_mean_result = train(cv_parameter_cell, matrix_cell_train, initial_matrix_cell)
%TRAIN Summary of this function goes here
%   Detailed explanation goes here
V1 = matrix_cell_train{1,1};
V2 = matrix_cell_train{2,1};
M  = matrix_cell_train{3,1};
%initial_matrix_cell = {W;H1;H2};
W = initial_matrix_cell{1,1};
H1 = initial_matrix_cell{2,1};
H2 = initial_matrix_cell{3,1};
%cv_parameter_cell = {alpha;beta; dimension_K_set; max_ites ; top_n_cv; cv_criteria};
alpha = cv_parameter_cell{1,1};
beta = cv_parameter_cell{2,1};
beta = cv_parameter_cell{2,1};
beta = cv_parameter_cell{2,1};
beta = cv_parameter_cell{2,1};
beta = cv_parameter_cell{2,1};
beta = cv_parameter_cell{2,1};
[L,W_out,H1_out,H2_out] =  CMNMF_LF(max_ites,InnerMaxIter,V1,V2,W,H1,H2,M,alpha,gama,lamta1,lamta2);
end

