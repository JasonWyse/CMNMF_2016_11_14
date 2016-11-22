function [] = Kernel(species_acronym,file_cell,parameter_cell)    
    test_groundTruth_str = parameter_cell{1,1};    
    g_p_network_file = file_cell{1,1};
    pathway_file_new = file_cell{2,1};
    ppi_file_new = file_cell{3,1};
    pathway_file_old = file_cell{4,1};
    ppi_file_old = file_cell{5,1};
    load(g_p_network_file,'g_p_network','gene_id');
    data_matrix = g_p_network;
    
    if strcmp(test_groundTruth_str,'ppi') == 1
      load(ppi_file_old,'ppi_old','ppi_old_index');
      [C,ia,ib] = intersect(gene_id,ppi_old_index);
      must_link_matrix = zeros(length(gene_id),length(gene_id));
      must_link_matrix(ia,ia) = ppi_old(ib,ib);
    elseif strcmp(test_groundTruth_str,'pathway') == 1
       load(pathway_file_old,'pathway_old');
       [groundTruth_pathway, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_old);
       groundTruth_pathway(groundTruth_pathway>1) = 1;
       [C,ia,ib] = intersect(gene_id,gene_idx_pathway);
       must_link_matrix = zeros(length(gene_id),length(gene_id));
       must_link_matrix(ia,ia) = groundTruth_pathway(ib,ib);
    end
    k = parameter_cell{2,1};
    method_dir = parameter_cell{3,1};
    centroid_matrix = zeros(k,size(data_matrix,2));
    cluster_matrix = zeros(size(data_matrix,1),k);
    [rows,columns] = find(must_link_matrix == 1);
     %1a
    count = 0;
    for i = 1:length(rows)
       row = rows(i);
       column = columns(i);
       if(row~=column)
          if(sum(cluster_matrix(row,:))==0&&sum(cluster_matrix(column,:))~=0)
              cluster_matrix(row,find(cluster_matrix(column,:)==1)) = 1;
          end
          if(sum(cluster_matrix(row,:))~=0&&sum(cluster_matrix(column,:))==0)
              cluster_matrix(column,find(cluster_matrix(row,:)==1)) = 1;
          end
          if(sum(cluster_matrix(row,:))==0&&sum(cluster_matrix(column,:))==0)
              new_neighb = find(sum(cluster_matrix)==0);
              %可能存在k给定值过小情况
              if count< k
                cluster_matrix(row,new_neighb(1)) =1 ;
                cluster_matrix(column,new_neighb(1)) = 1;
                count = count+1;
              end
          end
       end
    end
    
    lambda_neighb_index = find(sum(cluster_matrix)==0);
    if isempty(lambda_neighb_index)
        lambda_neighb = k;
    else
        lambda_neighb = lambda_neighb_index(1)-1;
    end
    %1b
    cluster_matrix_temp = cluster_matrix;
    [~,I] = sort(sum(cluster_matrix),2,'descend');
    for i = 1:k
       cluster_matrix(:,i) = cluster_matrix_temp(:,I(i)); 
    end
    %1c
    if lambda_neighb >= k
        for m = 1:k
            [rows,~] = find(cluster_matrix(:,m)==1);
            centroid_matrix(m,:) = mean(data_matrix(rows,:));
        end
    else
        for m = 1:lambda_neighb
            [rows,~] = find(cluster_matrix(:,m)==1);
            centroid_matrix(m,:) = mean(data_matrix(rows,:));
        end
        %把所有其他没有分类的点随机分类
        [rows,~] = find(sum(cluster_matrix,2)==0);
        left_neighb = k - lambda_neighb;
        avg_num = fix(length(rows)./left_neighb);
        for m = 1:left_neighb - 1
           cluster_matrix(rows(m*avg_num-avg_num+1:m*avg_num),m+lambda_neighb) = 1;
           [row,~] = find(cluster_matrix(:,m+lambda_neighb)==1);
           centroid_matrix(m+lambda_neighb,:) = mean(data_matrix(row,:));
        end
        cluster_matrix(rows((left_neighb-1)*avg_num+1:end),k) = 1;
        [row,~] = find(cluster_matrix(:,k)==1);
        centroid_matrix(k,:) = mean(data_matrix(row,:)); 
    end
    similarity_matrix_temp = (data_matrix * data_matrix');
    similarity_matrix = zeros(size(similarity_matrix_temp));
   for i = 1:size(similarity_matrix_temp,1)
        for j = 1:size(similarity_matrix_temp,1)
            similarity_matrix(i,j) = exp(-(similarity_matrix_temp(i,i) +similarity_matrix_temp(j,j)-2*similarity_matrix_temp(i,j))/18);
        end
    end
    K = similarity_matrix + must_link_matrix;
    
    Jpckm = calculate_J(K,cluster_matrix);
    step = 0;
    E = zeros(size(data_matrix,1),size(data_matrix,1));
    for i = 1:size(data_matrix,1)
        for j = 1:size(data_matrix,1)
            E(i,j) = similarity_matrix(i,i) +similarity_matrix(j,j)-2*similarity_matrix(i,j);
        end
    end
    while(1)
        tic;
        for i = 1:size(data_matrix,1)
           array = zeros(1,k);
           for j = 1:k
              [rows,~] = find(cluster_matrix(:,j) ==1);
              suma = 0;
              for m = 1:length(rows)
                  suma = suma + E(i,rows(m));
              end
              array(1,j) = suma/(length(rows));
           end
           [~,I] = max(array);
           cluster_matrix(i,:) = 0;
           cluster_matrix(i,I) = 1;
        end
        Jpckm2 = calculate_J(K,cluster_matrix);
        toc;
        if Jpckm2 >=Jpckm
            break;
        end
        Jpckm = Jpckm2;
        step = step+1;
    end
    
    [m,n]=size(g_p_network);
    F_measure_value = 1;
    G_predict = cluster_matrix;
    predict_GG_matrix = G_predict*G_predict';
    
    
    if strcmp(test_groundTruth_str,'pathway') ==1
        load(pathway_file_new, 'pathway_new');
        pathway_species_new =  pathway_new;
        [groundTruth_pathway, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_species_new);
        [F,Jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_id,groundTruth_pathway,...,
        gene_idx_pathway,F_measure_value);

        predicted_pathway_gene = transform(G_predict,gene_id);
        M_sim = get_Ratio(predicted_pathway_gene);
        Kernel_result = [F,Jaccard,RD,Precision,Recall,M_sim];
        Kernel_result_cell = {G_predict;Kernel_result};
        
        result_file_name = ['Kernel_' species_acronym '_result_pathway_' datestr(now,30) '.mat' ];  
        save([method_dir result_file_name],'Kernel_result_cell');
    end
    if strcmp(test_groundTruth_str,'ppi') ==1
       load(ppi_file_new, 'ppi_new','ppi_new_index');
        [F,Jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_id,ppi_new,...,
        ppi_new_index,F_measure_value);
    
        predicted_pathway_gene = transform(G_predict,gene_id);
        M_sim = get_Ratio(predicted_pathway_gene);
        Kernel_result = [F,Jaccard,RD,Precision,Recall,M_sim];
        Kernel_result_cell = {G_predict;Kernel_result};
        
        result_file_name = ['Kernel_' species_acronym '_result_ppi_' datestr(now,30) '.mat' ];  
        save([method_dir result_file_name],'Kernel_result_cell');
    end
    
    
end

function [Jpckm] = calculate_J(K,cluster_matrix)
    Jpckm = 0;
    Z = zeros(size(cluster_matrix,1),size(cluster_matrix,1));
    for i = 1:size(cluster_matrix,2)
       zc = cluster_matrix(:,i);
       Z(:,i) =  zc/sqrt(zc'*zc);
       Z(isnan(Z)) = 0;
    end
    Jpckm = sum(diag(Z'*K*Z));
end
function[predicted_pathway_gene] = transform(learned_matrix,gene_idx)
    predicted_pathway_gene = zeros(size(learned_matrix,2),max(sum(learned_matrix)));
    for i = 1:size(predicted_pathway_gene,1)
       [row,~] = find(learned_matrix(:,i)==1); 
       predicted_pathway_gene(i,1:length(row)) = gene_idx(row);
    end
    predicted_pathway_gene = predicted_pathway_gene';
end