function ia = getGeneSet(geneset, data, genes)

if(istable(data))
    if(~exist('genes','var') || isempty(genes))
        genes = data.Properties.RowNames;
    end
    data = data{:,:};
end

switch geneset   
    case 'table_s3'   
        % genes from Table S3
        table_s3 = load('data/table_s3.mat');
        selected_genes = table_s3.genes;
        [~,ia] = intersect(genes, selected_genes);
    case 'lm22'
        % LM22 genes
        lm22 = load('data/lm22.mat');
        selected_genes = lm22.genes;
        [~,ia] = intersect(genes, selected_genes);        
    case 'table_s12'   
        % genes from Table S12
        table_s12 = load('data/table_s12.mat');
        selected_genes = table_s12.genes;
        [~,ia] = intersect(genes, selected_genes);
    case 'all'
        % all genes
        ia = (1:size(data,1))';
    case 'merged'
         % merged Table S3, S12 and LM22 genes
        lm22 = load('data/lm22.mat');
        table_s3 = load('data/table_s3.mat');
        table_s12 = load('data/table_s12.mat');
        markers = load('data/marker_genes.mat');
        selected_genes = unique(union(lm22.genes, table_s3.genes));
        selected_genes = unique(union(selected_genes, table_s12.genes));
        selected_genes = unique(union(selected_genes, markers.genes));
        [~,ia] = intersect(genes, selected_genes);  
%     case 'tcell_genes'
% %         [~, sheets] = xlsfinfo('data/tg.xls');
% %         gene_list = {};
% %         for i=2:length(sheets)
% %             T = readtable('data/tg.xls','Sheet',sheets{i});
% %             gene_list = [gene_list; T{:,2}];
% %         end
% %         selected_genes = upper(unique(gene_list));
%         T = readtable('data/tcell_genes2.tsv','FileType','Text');
%         load('data/gene_info.mat')
%         tmp = innerjoin(T,gene_info,'LeftKeys','Gene',...
%             'RightKeys','EnsemblGeneID');
%         selected_genes = unique(tmp.ApprovedSymbol);
%         table_s12 = load('data/table_s12.mat');
%         selected_genes = unique(union(table_s12.genes, selected_genes));
%         [~,ia] = intersect(genes, selected_genes);
%     case 'merged'
%         % merged Table S3 and LM22 genes
%         lm22 = load('data/lm22.mat');
%         table_s3 = load('data/table_s3.mat');
%         selected_genes = unique(union(lm22.genes, table_s3.genes));
%         [~,ia] = intersect(genes, selected_genes);        
%     case 'cutpredictors'
%         load('cutpredictor.mat');
%         selected_genes = strrep(cutpredictors(2:end),'_','-');
%         [~,ia] = intersect(genes, selected_genes); 
%     case 'cutpredictors_pca'
%         load('data/cutpredictors_pca.mat');
%         selected_genes = cutpredictors;
%         [~,ia] = intersect(genes, selected_genes); 
end

% keep only genes that are expressed in at least one column
ib = any(data(ia,:)>0,2);
ia = ia(ib);
