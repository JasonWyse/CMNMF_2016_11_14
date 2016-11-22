function [learned_matrix_cell] =  HNMF_Train(train_parameter_cell, matrix_cell_train, initial_matrix_cell)
    V1 = matrix_cell_train{1,1};
    V2 = matrix_cell_train{2,1};
    M  = matrix_cell_train{3,1};
    %initial_matrix_cell = {W;H1;H2};
    W = initial_matrix_cell{1,1};
    H1 = initial_matrix_cell{2,1};
    H2 = initial_matrix_cell{3,1};
    % cv_parameter_cell = {alpha; beta; max_ites; top_n_cv; cv_criteria};
    lambda1 = train_parameter_cell{1,1};
    lambda2 = train_parameter_cell{2,1};
    MaxIter = train_parameter_cell{3,1};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %[W,H,tempV,rmse,isCoverage] = HNMF_sparse(Max_iter,V1,V2,W,H1,H2,M,lambda1,lambda2)
    [W_out,H_out,V_out,~,isCoverage] = HNMF_sparse(MaxIter,V1,V2,W,H1,H2,M,lambda1,lambda2);
    learned_matrix_cell = {isCoverage; W_out; H_out; V_out};

end

