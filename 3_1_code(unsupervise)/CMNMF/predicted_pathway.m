function [Z_filtered,predicted_pathway_gene] = predicted_pathway(W,K,mgi_id)
%计算Z-score
z_threshold =K;
Z=zscore(W');
Z(Z<z_threshold)=0;
Z =Z';

%聚簇
[m,n]=size(Z);
%k(每个基因，每次筛选的个数）
k = zeros(m,1);
%Z_filtered为预测出的gene-pathway关系矩阵，为0,1矩阵
Z_filtered=zeros(m,n);
Z(Z>0) = 1;
Z_filtered = Z;
    
%predicted_pathway_mgi_id每一列为一个pathway,其中包含此pathway中基因的mgi_id
%将Z_filter矩阵转化为predicted_pathway_mgi_id矩阵

max_pathway_gene = max(sum(Z_filtered));
gene_location = zeros(max_pathway_gene,n);
predicted_pathway_gene = zeros(max_pathway_gene,n);
     
   for i = 1:n
         gene_location = find(Z_filtered(:,i)==1);
         gene_mgi_id = mgi_id(gene_location);
         predicted_pathway_gene(1:length(gene_mgi_id),i)=gene_mgi_id;
   end
   

     
