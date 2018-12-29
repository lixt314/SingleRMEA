function [celltype_matrix, celltype_expression, cellnames, celltype] = celltype_markers(data,tsneX,id_celltype,thresh,id)

genes = {'Unknown',{},{},{};...
    'T cells',{'CD3D','CD3E','CD3G','CD27','CD28'},{},{};...
    'CD4+ T cells',{'CD4'},{},{'FOXP3','IL2RA','CTLA4'};...
    'CD8+ T cells',{'CD8B'},{'CD8A'},{'CD4'};...
    'regulatory T cells',{'FOXP3','IL2RA','CD4','CTLA4'},{},{};...
    'B cells',{'CD19','MS4A1','CD79A','CD79B','BLNK'},{},{};...
    'Macrophages/Monocytes',{'CD14','CD68','CD163','CSF1R','FCGR3A'},{},{};...
    'Dendritic cells',{'IL3RA','CLEC4C','NRP1'},{},{};...
    'Natural killer cells',{'FCGR3A','FCGR3B','NCAM1','KLRB1','KLRB1','KLRC1','KLRD1','KLRF1','KLRK1'},{},{};...
    'Endothelial cells',{'VWF','CDH5','SELE'},{},{};...
    'Cancer associated fibroblasts',{'FAP','THY1','COL1A1','COL3A1'},{},{};...
    'Ovarian carcinoma cells',{'WFDC2','EPCAM','MCAM'},{},{};...
    'Melanoma cells',{'PMEL','MLANA','TYR','MITF'},{},{}}; % J Clin Aesthet Dermatol. 2014 Jun; 7(6): 13?24.

if(~exist('id_celltype','var') || isempty(id_celltype))
    id_celltype = 1:size(genes,1);
end
if(~exist('thresh','var') || isempty(thresh))
    thresh = 0;
end
if(~exist('id','var') || isempty(id))
    id = true(1,size(tsneX,1));
end

lb = min(tsneX);
ub = max(tsneX);

%% reduce to id
data = data(:,id);
tsneX = tsneX(id,:);

%%
celltype_matrix = zeros(size(genes,1),size(data,2));
celltype_expression = zeros(size(genes,1),size(data,2));
celltype = zeros(1,size(data,2));
cellnames = genes(:,1)';

if length(id_celltype)<10
    m = 3;
else
    m = 4;
end

for i=1:length(id_celltype)
    j = id_celltype(i);
    
    [~, idx] = intersect(data.Properties.RowNames,genes{j,2});
    or_marker_expr = data{idx,:};
    or_marker_expr = bsxfun(@rdivide,...
        or_marker_expr,...
        max(or_marker_expr,[],2)); 
    or_marker_expr = mean(or_marker_expr,1);
    
    if ~isempty(genes{j,3})
        [~, idx] = intersect(data.Properties.RowNames,genes{j,3});
        and_marker_expr = data{idx,:};
        and_marker_expr = bsxfun(@rdivide,...
            and_marker_expr,...
            max(and_marker_expr,[],2));
        and_marker_expr = mean(and_marker_expr,1);
    else
        and_marker_expr = ones(size(or_marker_expr));
    end
    
    if ~isempty(genes{j,4})
        [~, idx] = intersect(data.Properties.RowNames,genes{j,4});
        not_marker_expr = data{idx,:};
        not_marker_expr = bsxfun(@rdivide,...
            not_marker_expr,...
            max(not_marker_expr,[],2));

        not_marker_expr = 1-mean(not_marker_expr,1);
    else
        not_marker_expr = ones(size(or_marker_expr));
    end
    
    marker_expr = or_marker_expr.*and_marker_expr.*not_marker_expr;
    
    marker_expr = bsxfun(@rdivide,marker_expr,max(marker_expr,[],2));
    
    celltype_matrix(j,:) = marker_expr>thresh;
    celltype_expression(j,:) = marker_expr;
    celltype(marker_expr>thresh) = j-1;
    
    [marker_expr, id_sort] = sort(marker_expr);
    
    subaxis(m,m,i,'S',0.02,'SV',0.03,'M',0.03,'MT',0.05)
    s = scatter(tsneX(id_sort,1),tsneX(id_sort,2),5,marker_expr);
    s.MarkerFaceColor = s.MarkerEdgeColor;
    s.MarkerFaceAlpha = 0.5;
    s.MarkerEdgeAlpha = 0.5;
    ax = gca;
    colormap grayorangered;
    ax.XTick = [];
    ax.YTick = [];
    xlim([-1.05*abs(lb(1)) 1.05*abs(ub(1))]);
    ylim([-1.05*abs(lb(2)) 1.05*abs(ub(2))]);
    box on;
    title(genes{j,1})    
end
