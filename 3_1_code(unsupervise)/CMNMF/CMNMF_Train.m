function [learned_matrix_cell] =  CMNMF_Train(train_parameter_cell, matrix_cell_train, initial_matrix_cell)
    V1 = matrix_cell_train{1,1};
    V2 = matrix_cell_train{2,1};
    M  = matrix_cell_train{3,1};
    %initial_matrix_cell = {W;H1;H2};
    W = initial_matrix_cell{1,1};
    H1 = initial_matrix_cell{2,1};
    H2 = initial_matrix_cell{3,1};
    % cv_parameter_cell = {alpha; beta; max_ites; top_n_cv; cv_criteria};
    alpha = train_parameter_cell{1,1};
    beta = train_parameter_cell{2,1};
    max_ites = train_parameter_cell{3,1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    lambda1 = 0;
    lambda2 = 0;
    [L,W_out,H1_out,H2_out] =  CMNMF_LF(max_ites,V1,V2,W,H1,H2,M,alpha,beta,lambda1,lambda2);
    learned_matrix_cell = {L; W_out; H1_out; H2_out};

end

