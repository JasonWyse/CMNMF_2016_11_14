function [ F,jaccard,RD,Precision,Recall ] = rand_index( predict,predict_row,standard,standard_row,a)

[~,ia,ib] = intersect(predict_row,standard_row);
g_g_predict = predict(ia,ia);
g_g_standard = standard(ib,ib);

g_g_predict(g_g_predict>0) = 1;
g_g_predict = g_g_predict-diag(ones(size(g_g_predict,1),1));
g_g_standard(g_g_standard>0) = 1;
g_g_standard = g_g_standard-diag(ones(size(g_g_standard,1),1));
g_g_predict(g_g_predict<0) = 0;
g_g_standard(g_g_standard<0) = 0;

TP = nnz(g_g_predict+g_g_standard==2);
TN = nnz(g_g_predict+g_g_standard==0);
FP = nnz(g_g_predict-g_g_standard==1);
FN = nnz(g_g_predict-g_g_standard==-1);

RD = (TP+TN)/(TP+TN+FP+FN);

Precision = TP/(TP+FP);

Recall = TP/(TP+FN);

F=(a^2+1)*Recall*Precision/(a^2*(Recall+Precision));

jaccard = TP/nnz(g_g_predict+g_g_standard);

end

