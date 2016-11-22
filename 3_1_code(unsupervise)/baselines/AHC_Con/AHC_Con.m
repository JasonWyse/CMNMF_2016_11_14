function [] = AHC_Con(species_acronym,file_cell,parameter_cell)
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
    kend = parameter_cell{2,1};
    method_dir = parameter_cell{3,1};
    k=400;
    centroid_matrix = zeros(k,size(data_matrix,2));
    cluster_matrix = zeros(size(data_matrix,1),k);
    [rows,columns] = find(must_link_matrix == 1);
     %1a    生成must_link_matrix
    for i = 1:k
       cluster_matrix(i,i) = 1; 
       centroid_matrix(i,:) = data_matrix(i,:);
    end
    for i = 1:size(must_link_matrix,1)
        must_link_matrix(i,i) = 1; 
    end
    must_link_matrix = must_link_matrix + must_link_matrix';
    must_link_matrix(must_link_matrix>1) = 1;
    step = 1;
    % 先把前面400个点当做400个簇，然后将后面的点一个个放入400簇中
    while step >0
       step = step -1;
       tic;
       for i = k+1:size(data_matrix,1)
          %计算该点距离400簇哪一簇最近
          self_matrix = zeros(k,size(data_matrix,2));
          for j = 1:k
             self_matrix(j,:) = data_matrix(i,:); 
          end
          array = sum((self_matrix-centroid_matrix).^2,2);
          %得到距离矩阵后，根据must_link约束找又满足约束又最近的一簇
          [~,columns] = find(must_link_matrix(i,:)==1);
          [~,clusters_good] = find(cluster_matrix(columns,:)==1);
          [~,I] = min(array(clusters_good));
          I = clusters_good(I);
          if isempty(I)==1
              [~,I] = min(array);
          end
          I = I(1);
          cluster_matrix(i,:) = 0;
          cluster_matrix(i,I) = 1;
       end
       for i = 1:k
          [rows,~] = find(cluster_matrix(:,i)==1);
          if(length(rows)==1)
              centroid_matrix(i,:) = data_matrix(rows,:);
          else
              centroid_matrix(i,:) = sum(data_matrix(rows,:))/length(rows);
          end
       end   
       toc;
    end
    
    %算出一开始的簇差矩阵
    dissimiliarty_matrix = size(k,k);
    for i = 1:k
        for j = i:k
            result = sum((centroid_matrix(i,:)-centroid_matrix(j,:)).^2);
            dissimiliarty_matrix(i,j) = result;
            dissimiliarty_matrix(j,i) = result;
        end
    end
    while kend < k
        for i = 1:k
            dissimiliarty_matrix(i,i) = inf;
        end
        %找到现在的k簇中最接近的两个簇
        m = min(min(dissimiliarty_matrix));
        [rows,columns] = find(dissimiliarty_matrix ==m);
        Group1 = rows(1);
        Group2 = columns(1);
        if Group1>Group2
            Group1 = columns(1);
            Group2 = rows(1);
        end
        %把k*k矩阵变为k-1*k-1矩阵
        dissimiliarty_matrix_temp = zeros(k-1,k-1);
        dissimiliarty_matrix_used = dissimiliarty_matrix;
        dissimiliarty_matrix_used([Group1,Group2],:) = [];
        dissimiliarty_matrix_used(:,[Group1,Group2]) = [];
        dissimiliarty_matrix_temp(1:end-1,1:end-1) = dissimiliarty_matrix_used;
        
        [rows1,~] = find(cluster_matrix(:,Group1)==1);
        [rows2,~] = find(cluster_matrix(:,Group2)==1);
        Gp = length(rows1);
        Gq = length(rows2);
        cluster_matrix(:,[Group1,Group2]) =[];
        cluster_matrix(:,k-1) = 0;
        cluster_matrix([rows1;rows2],k-1) =1;
        for i = 1:k-2
            [rows3,~] = find(cluster_matrix(:,i)==1);
            G = length(rows3);
            if i < Group1
                distance = Gp/(Gp+Gq)*dissimiliarty_matrix(Group1,i)+...,
                Gq/(Gp+Gq)*dissimiliarty_matrix(Group2,i)-...,
                Gp*Gq/((Gp+Gq).^2)*dissimiliarty_matrix(Group1,Group2);
            elseif i>Group1 && i<Group2
                distance = Gp/(Gp+Gq)*dissimiliarty_matrix(Group1,i+1)+...,
                Gq/(Gp+Gq)*dissimiliarty_matrix(Group2,i+1)-...,
                Gp*Gq/((Gp+Gq).^2)*dissimiliarty_matrix(Group1,Group2);
            elseif i >Group2
                distance = Gp/(Gp+Gq) *dissimiliarty_matrix(Group1,i+2)+...,
                Gq/(Gp+Gq)*dissimiliarty_matrix(Group2,i+2)-...,
                Gp*Gq/((Gp+Gq).^2)*dissimiliarty_matrix(Group1,Group2);
            end
            dissimiliarty_matrix_temp(k-1,i) = distance;
            dissimiliarty_matrix_temp(i,k-1) = distance;
        end
        dissimiliarty_matrix = dissimiliarty_matrix_temp;
        k = k-1;
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
        AHC_Con_result = [F,Jaccard,RD,Precision,Recall,M_sim];
        AHC_Con_result_cell = {G_predict;AHC_Con_result};
        
        result_file_name = ['AHC_Con_filled_' species_acronym '_result_pathway_' datestr(now,30) '.mat' ];  
        save([method_dir result_file_name],'AHC_Con_result_cell');
    end
    if strcmp(test_groundTruth_str,'ppi') ==1
        %load(ppi_file_new, 'ppi_new','ppi_new_index');
        load(ppi_file_new, 'ppi_difference','ppi_new_index');
        ppi_new = ppi_difference;
        [F,Jaccard,RD,Precision,Recall] = rand_index( predict_GG_matrix,gene_id,ppi_new,...,
        ppi_new_index,F_measure_value);
    
        predicted_pathway_gene = transform(G_predict,gene_id);
        M_sim = get_Ratio(predicted_pathway_gene);
        AHC_Con_result = [F,Jaccard,RD,Precision,Recall,M_sim];
        AHC_Con_result_cell = {G_predict;AHC_Con_result};
        
        result_file_name = ['AHC_Con_filled_' species_acronym '_result_ppi_difference_' datestr(now,30) '.mat' ];  
        save([method_dir result_file_name],'AHC_Con_result_cell');
    end
end
function[predicted_pathway_gene] = transform(learned_matrix,gene_id)
    predicted_pathway_gene = zeros(size(learned_matrix,2),max(sum(learned_matrix)));
    for i = 1:size(predicted_pathway_gene,1)
       [row,~] = find(learned_matrix(:,i)==1); 
       predicted_pathway_gene(i,1:length(row)) = gene_id(row);
    end
    predicted_pathway_gene = predicted_pathway_gene';
end