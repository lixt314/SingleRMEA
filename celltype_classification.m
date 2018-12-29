function [training_celltype]=celltype_classification(pbmc_dataset,geneset,savefigure)
addpath('tools')
addpath('tools/bhtsne/')


%%
fname = sprintf('data/normalized_data_%s.mat', pbmc_dataset)
load(fname);


fname = sprintf('data/classified_data_%s_%s.mat', pbmc_dataset, geneset)

if 1%(~exist(fname,'file'))
    %% restrict to 'merged' gene set
    ia = getGeneSet(geneset,data);
 

    %% run tSNE
    fname2 = sprintf('data/tsne_data_%s_%s.mat', pbmc_dataset, geneset);
    if (~exist(fname2,'file'))
         addpath('tools/bhtsne/');
        tsneX = fast_tsne(data{ia,:}');
        save(fname2,'tsneX');
    else
        load(fname2);
    end
   
    %% DBSCAN
    if strcmp(pbmc_dataset,'4k')
        epsilon = 2;
        ncells = 25;
    elseif strcmp(pbmc_dataset,'8k')
        epsilon = 2;
        ncells = 25;
    elseif strcmp(pbmc_dataset,'12k')
        epsilon = 2;
        ncells = 30;
    elseif strcmp(pbmc_dataset,'donor_abc')
        epsilon = 0.7;
        ncells = 30;
    elseif strcmp(pbmc_dataset,'donor_abc_12k')
        epsilon = 1.5;
        ncells = 25;
    end
   %%%%%%DBSCAN
    
  %  cluster = DBSCAN(tsneX,epsilon,ncells);

  %  cluster(cluster==-1) = 0;
    %CL=max(cluster)
  % silhouette(tsneX, cluster)

      dist=pdist2(tsneX,tsneX);
   
    [cluster] = cluster_dp(dist);
    cluster=cluster';
   % cluster=cluster-1;
       cluster(cluster==-1) = 0; 
   %    cluster1=cluster';
   
   % XX = matpy.mat2nparray(tsneX);
   % labels= matpy.mat2nparray(cluster);
   % SC=py.sklearn.metrics.silhouette_score(XX,labels)
    %CH=py.sklearn.metrics.calinski_harabaz_score(XX,labels)
  

    
    load('colors.mat');

  
    dbscan_cols = parula(length(unique(cluster)));
    rng(0);
    dbscan_cols = dbscan_cols(randperm(max(unique(cluster))),:);
    dbscan_cols = [0.8*ones(1,3); dbscan_cols];

    fig = figure(1);
    clf
    lims = [min(tsneX(:)) max(tsneX(:))];
    

    for ci=1:length(indication)
        id = ~cellfun(@isempty,strfind(sample_id,indication{ci}));
  
        my_gscatter(tsneX(id,1),tsneX(id,2),cluster(id),dbscan_cols,symbs(ci),false,[],25,lims);
    end

    my_gscatter(tsneX(:,1),tsneX(:,2),cluster,[],{'none'},true,[],25,lims);

    fig.Color = 'w';
    fig.Position = [0 0 900 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
  
    fname = sprintf('figures/dbscan_clusters_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        %     fig.Renderer = 'painters';
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
        %     print(fig,'-dpdf',[fname '.pdf'],'-r300');
    end

    %% marker genes
    fig = figure(2);
    clf
    [celltype_matrix, celltype_expression, cellnames] = ...
        celltype_markers(data,tsneX,[2 6:13],0.3); 
    fig.Color = 'w';
    fig.Position = [0 0 800 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
    fname = sprintf('figures/marker_expression_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        %     fig.Renderer = 'painters';
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
        %     print(fig,'-dpdf',[fname '.pdf'],'-r300');
    end

    %% training data
    training_celltype = assign_clusters(cluster,celltype_matrix,celltype_expression);

    % remove bad training data
    id_endo = 9; id_caf = 10; id_ov = 11; id_mel = 12; id_nk = 8;
    id_macro = 6; id_bcell = 5; id_unknown = 0;
    id_melanoma = findExpr(sample_id,'melanoma');
    id_ascites = findExpr(sample_id,'ascites');
    id_pbmc = findExpr(sample_id,'pbmc');
    training_celltype(training_celltype==id_ov & id_melanoma) = id_unknown;
    training_celltype(training_celltype==id_mel & id_ascites) = id_unknown;
    training_celltype(training_celltype==id_endo & id_pbmc) = id_unknown;
    training_celltype(training_celltype==id_caf & id_pbmc) = id_unknown;
    training_celltype(training_celltype==id_ov & id_pbmc) = id_unknown;
    training_celltype(training_celltype==id_mel & id_pbmc) = id_unknown;
    fprintf(1,'%d/%d celltypes found for classification training\n',sum(unique(training_celltype)~=0),9);
    id = setdiff([0:1 5:12],unique(training_celltype))+1

    if ~isempty(id)
        cellnames(id)
        return;
    end
   %save('data/training_celltype.mat','training_celltype')
 
    %% plot training set
    fig = figure(3);
    clf;

    lims = [min(tsneX(:)) max(tsneX(:))];
    for ci=1:length(indication)
        id = ~cellfun(@isempty,strfind(sample_id,indication{ci}));
        my_gscatter(tsneX(id,1),tsneX(id,2),training_celltype(id),cols,symbs(ci),false,[],25,lims);
    end

    fig.Color = 'w';
    fig.Position = [0 0 900 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
    fname = sprintf('figures/celltype_training_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        %     fig.Renderer = 'painters';
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
        %     print(fig,'-dpdf',[fname '.pdf'],'-r300');
    end

    %% classify
%size(training_celltype)

  % [predicted_celltype, trainedClassifier] = classify_celltypes(data(ia,:),training_celltype);
 [predicted_celltype] = Ensemble_Classify(data{ia,:}',training_celltype);

    % remove bad classifications
    id_endo = 9; id_caf = 10; id_ov = 11; id_mel = 12; id_nk = 8;
    id_macro = 6; id_bcell = 5; id_unknown = 0;
    id_melanoma = findExpr(sample_id,'melanoma');
    id_ascites = findExpr(sample_id,'ascites');
    id_pbmc = findExpr(sample_id,'pbmc');
    predicted_celltype(predicted_celltype==id_ov & id_melanoma) = id_unknown;
    predicted_celltype(predicted_celltype==id_mel & id_ascites) = id_unknown;
    predicted_celltype(predicted_celltype==id_endo & id_pbmc) = id_unknown;
    predicted_celltype(predicted_celltype==id_caf & id_pbmc) = id_unknown;
    predicted_celltype(predicted_celltype==id_ov & id_pbmc) = id_unknown;
    predicted_celltype(predicted_celltype==id_mel & id_pbmc) = id_unknown;

    
    %size(find(predicted_celltype==training_celltype),2)./size(predicted_celltype,2)
   % fas
    %% plot classification result
    fig = figure(4);
    clf;
    lims = [min(tsneX(:)) max(tsneX(:))];
    for ci=1:length(indication)
        id = ~cellfun(@isempty,strfind(sample_id,indication{ci}));
        my_gscatter(tsneX(id,1),tsneX(id,2),predicted_celltype(id),cols,symbs(ci),false,[],25,lims);
    end

    fig.Color = 'w';
    fig.Position = [0 0 900 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
    fname = sprintf('figures/celltype_classified_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        %     fig.Renderer = 'painters';
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
        %     print(fig,'-dpdf',[fname '.pdf'],'-r300');
    end

    %%

end
