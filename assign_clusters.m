function training_celltype = assign_clusters(cluster,celltype_matrix,celltype_expression)

% calculate score for each cluster and select highest ranked celltype
unique_cluster = unique(cluster);

cluster_score = NaN(size(celltype_matrix,1),length(unique_cluster));
cluster_celltype = cell(1,size(celltype_matrix,1));
for i=1:length(unique_cluster)
    cluster_score(:,i) = sum(celltype_expression(:,cluster==unique_cluster(i)),2);
    ncellspercluster = sum(celltype_matrix(:,cluster==unique_cluster(i)),2);
    cluster_score(ncellspercluster<5,i) = 0;
    [sorted_score, id] = sort(cluster_score(:,i),'descend');
    % select top 25% scores
    id = id(sorted_score>0.75*max(sorted_score));
    for j=1:length(id)
        cluster_celltype{id(j)} = [cluster_celltype{id(j)}, unique_cluster(i)];
    end
end

% select cells in each assigned cluster
training_celltype = zeros(1,length(cluster));
for i=1:size(celltype_matrix,1);
    cluster_id = cluster_celltype{i};
    id = false(1,length(cluster));
    for j=1:length(cluster_id)
        id = id | (celltype_matrix(i,:) & cluster'==cluster_id(j));
    end
    training_celltype(id) = i-1;
end
