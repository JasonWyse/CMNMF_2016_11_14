function [learned_matrix_cell] = NMF_Train(cv_parameter_cell, matrix_cell_train, initial_matrix_cell)%( MaxIter,V,W,H)   
    %L为目标函数值
    %initial_matrix_cell = {W;H1;H2};matrix_cell_train = {V1; V2; M};
    MaxIter = cv_parameter_cell{3,1};
    W = initial_matrix_cell{1,1};
    H1 = initial_matrix_cell{2,1};
    H2 = initial_matrix_cell{3,1};
    H = [H1,H2];
    V1 = matrix_cell_train{1,1};
    V2 = matrix_cell_train{2,1};
    M = matrix_cell_train{3,1};
    first_level_id = matrix_cell_train{4,1};
    V = [V1,V2];
    V = Extend_GPNetwork_TruePathRule( V, M, first_level_id);
    L=zeros(MaxIter*2,1);
     L(1,1)=objective(V,W,H);
    for i=1:MaxIter
       %update W             
        W=W.*(V*H')./(W*(H*H')+eps);
        t = (i-1)*2+2;
        L(t,1)=objective(V,W,H);
        if(L(t,1) > L(t-1,1))
          break;
        end
        W(isnan(W))=0; 
        %update H    
        H=H.*(W'*V)./(W'*W*H+eps);
        t = (i-1)*2+3;
        L(t,1)=objective(V,W,H);
        if L(t,1)>L((t-1),1)
            break;
        end            
        H(isnan(H))=0;    
    end       
    
    D = diag(1./sqrt(sum(W.^2,1)));%W矩阵每个元素平方后，每列的和构成的向量作为对角矩阵对角线上元素
    D_inverse = diag(sqrt(sum(W.^2,1)));
    W = W*D;
    H = D_inverse * H;    
    learned_matrix_cell = {L;W;H};
end

function [O]= objective(V,W,H)
    
    O1=sum(sum((V-W*H).^2));
    %O2=-gama*trace(H*M*H2')+lamat1*sum(sum((W).^2))+lamat2*(sum(sum(H.^2))+sum(sum(H2.^2)));
    O=O1;    
end


