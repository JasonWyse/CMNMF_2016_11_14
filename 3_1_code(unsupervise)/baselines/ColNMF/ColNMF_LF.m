function [L,W,H1,H2] = ColNMF_LF( MaxIter,A1,A2,W,H1,H2,M,alpha,beta,lambda1,lambda2)
%L为目标函数值
    D1 = diag(sum(M,2));
    D2 = diag(sum(M,1));
 %L=zeros(1+2*Inner_MaxIter,MaxIter);
 L=zeros(MaxIter*3,3);
 [L(1,1),L(1,2),L(1,3)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
for i=1:MaxIter  
    %update W 
    W=W.*(A1*H1'+alpha*A2*H2')./(W*(H1*H1')+alpha*W*(H2*H2')+lambda1*W+eps);
    t = (i-1)*3+2;
    [L(t,1),L(t,2),L(t,3)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    if(L(t,1) > L(t-1,1))
      break;
    end
    W(isnan(W))=0;    
    %update H1
    H1=H1.*((W'*A1 + beta*H2*M')./(W'*W*H1 + beta*H1*D1 + eps));
    t = (i-1)*3+3;
    [L(t,1),L(t,2),L(t,3)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    if L(t,1)>L((t-1),1)
        break;
    end            
    H1(isnan(H1))=0; 
    %update H2         
    H2=H2.*((alpha*W'*A2 + beta*H1*M)./(alpha*(W'*W*H2) + beta*H2*D2 + eps));  
    t = (i-1)*3+4;
    [L(t,1),L(t,2),L(t,3)]=objective(alpha,beta,lambda1,lambda2,A1,A2,W,H1,H2,M);
    if L(t,1)>L((t-1),1)
        break;
    end          
    H2(isnan(H2))=0; 
end
D = diag(1./sqrt(sum(W.^2,1)));%W矩阵每个元素平方后，每列的和构成的向量作为对角矩阵对角线上元素
D_inverse = diag(sqrt(sum(W.^2,1)));
W = W*D;
H1 = D_inverse * H1;
H2 = D_inverse * H2;
% W=diag(1./sum(W,2))*W;
% W(isnan(W))=0;
% H1=H1*diag(1./sum(H1));
% H2=H2*diag(1./sum(H2));
end

function [O,O1,O2]= objective(alpha,beta,lambda1,lambda2,V1,V2,W,H1,H2,M)
    
    O1=sum(sum((V1-W*H1).^2))+alpha*sum(sum((V2-W*H2).^2));
    %O2=-gama*trace(H1*M*H2')+lamat1*sum(sum((W).^2))+lamat2*(sum(sum(H1.^2))+sum(sum(H2.^2)));
    D1 = diag(sum(M,2));
    D2 = diag(sum(M,1));
    O2=beta*(trace(H1*D1*H1') + trace(H2*D2*H2')-2*trace(H1*M*H2'));
    O=O1+O2;
    
end


