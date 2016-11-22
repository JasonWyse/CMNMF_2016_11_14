function [learned_matrix_cell,best_parameter_array,evaluation_result] = learn(input_parameter_cell, matrix_cell_train,...,
    matrix_cell_validation, initial_matrixFileName_cell)
    alpha_set = input_parameter_cell{1,1};
    beta_set = input_parameter_cell{2,1};
    max_ites = input_parameter_cell{3,1};
    top_n_set = input_parameter_cell{4,1};
    cv_criteria = input_parameter_cell{5,1};
    method_dir = input_parameter_cell{6,1};
    species = input_parameter_cell{7,1};
    top_n_cv = top_n_set(end);
    i=1;
    %each row:[RD,F,Precision,Recall,jaccard]
    %evaluation_result:每组参数组合下，十次随机初始化学习到十次结果的平均值
    evaluation_result = zeros(length(alpha_set)*length(beta_set),8);

    for alpha = alpha_set
        for beta = beta_set
            cv_parameter_cell = {alpha; beta; max_ites; top_n_cv; cv_criteria; method_dir; species};
            [evaluation_result(i,1:end-2),~] = cv_train(cv_parameter_cell, matrix_cell_train, matrix_cell_validation,...,
                initial_matrixFileName_cell);
            evaluation_result(i,end-1:end) = [alpha,beta];
            i = i+1;
        end
    end

    best_parameter_array = get_best_parameter(evaluation_result,cv_criteria);
    %获取到最优参数后，我们选择十次随机初始化结果中，性能最好的一次初始化结果最为最后评价用的结果
    cv_parameter_cell = {best_parameter_array(1); best_parameter_array(2); max_ites; top_n_cv; cv_criteria; method_dir; species};
    [~,evaluation_result_best_parameter] = cv_train(cv_parameter_cell, matrix_cell_train, matrix_cell_validation,...,
                initial_matrixFileName_cell);
    %evaluation_result = [RD,F,Precision,Recall,jaccard],F1 value is in the
    %second column
    [M,I]  = max(evaluation_result_best_parameter(:,1));
    method_data_dir = [method_dir 'data/' species '/'];
    file_name = [method_data_dir initial_matrixFileName_cell{I,1}];
    %file_name = initial_matrixFileName_cell{I,1};
    load(file_name);
    initial_matrix_cell = {W;H1;H2};
    [learned_matrix_cell] =  NMF_Train(cv_parameter_cell, matrix_cell_train, initial_matrix_cell);
end
%cv_train 返回在一组参数组合情况下，十次随机初始化学习到的十次结果
function  [evaluation_result_average,evaluation_result] = cv_train(cv_parameter_cell, matrix_cell_train, matrix_cell_validation, initial_matrixFileName_cell)
    train_parameter_cell = cv_parameter_cell;
    method_dir = cv_parameter_cell{6,1};
    species = cv_parameter_cell{7,1};
    method_data_dir = [method_dir 'data/' species '/'];
    %given a fixed parameters, return the average result of ten times' different initial matrix values
    %initial_matrix_cell = {file_name_cell;W; H1; H2};
    evaluation_result = zeros(length(initial_matrixFileName_cell),6);
    for i=1:length(initial_matrixFileName_cell)
        file_name = [method_data_dir initial_matrixFileName_cell{i,1}];
        load(file_name);
        initial_matrix_cell = {W;H1;H2};
        [learned_matrix_cell] = NMF_Train(train_parameter_cell, matrix_cell_train, initial_matrix_cell);    
        evaluation_result(i,:) = NMF_Evaluate(learned_matrix_cell, matrix_cell_validation,train_parameter_cell);
        
    end
    evaluation_result_average = mean(evaluation_result,1);
end
function[best_parameter_array] = get_best_parameter(evaluation_result,cv_criteria)
% each row:[F,Jaccard,RD,Precision,Recall] [RD,F,Precision,Recall,Jaccard]
    if strcmp(cv_criteria,'F1') == 1
     [M,I]  = max(evaluation_result(:,1));
     best_parameter_array = evaluation_result(I,end-1:end);
    elseif strcmp(cv_criteria,'Jaccard') == 1
     [M,I]  = max(evaluation_result(:,2));
     best_parameter_array = evaluation_result(I,end-1:end);
    elseif strcmp(cv_criteria,'RD') == 1
     [M,I]  = max(evaluation_result(:,3));
     best_parameter_array = evaluation_result(I,end-1:end);   
    elseif strcmp(cv_criteria,'Precision') == 1
     [M,I]  = max(evaluation_result(:,4));
     best_parameter_array = evaluation_result(I,end-1:end); 
     elseif strcmp(cv_criteria,'Recall') == 1
     [M,I]  = max(evaluation_result(:,5));
     best_parameter_array = evaluation_result(I,end-1:end);
    end

end