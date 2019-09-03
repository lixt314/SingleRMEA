# SingleRMEA: Elucidating Gene Expression of Cells from Single-Cell RNA-seq data

The development of single-cell RNA-seq technologies provides new opportunities for biology since it has become an accepted experimental method throughput improvements enabling applications for cell type discovery. However,  high-throughput applications of single-cell RNA-seq to solid tissues rely on the formal cell type definitions. Unfortunately, it is unclear how to formulate such definitions to analyze the crucial information on cells' location further since high levels of technical noise in most data. To address this challenge, we present a computational method (SingleRMEA) to conduct the robust cell type classifiers including clustering and classification to address large-scale PBMC dataset and human tissue sources composed of complex mixtures of cell types and subtypes. For clustering, clustering by fast search and find of density peaks (CDP) is employed to perform clustering to partition the cells into a few distinct subpopulations. For classification, we develop a new ensemble construction method to predict the cell type for single-cell RNA-seq, which applies multiobjective optimization to the stacking ensemble construction process to generate domain-specific configurations under two objective functions. To validate our SingleRMEA method, we compare its performance across two PBMCs datasets including PBMC-4k and PBMC-12 merged data. The experimental results demonstrate that the SingleRMEA can obtain superior performance over the current state-of-the-art methods. Meanwhile, it also demonstrates that SingleRMEA can enable the construction of cell type classifiers that can be direct to other new single-cell RNA-seq data. SingleRMEA is written in Matlab and available at https://github.com/lixt314/SingleRMEA.

# Data

We have collected the immune cell single-cell RNA-seq profiles across the human datasets from \citet{schelker2017estimation}. In particular, two PBMC-merged data sources are considered by interpreting three different PBMCs data with the melanoma patient samples and ovarian cancer ascites samples. The PBMCs data includes PBMC-4k and PBMC-12k. The first dataset (PBMC-4k merged) includes 4k PBMCs \footnote{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k} from a healthy donor, which contains 4000 single cells sequenced on Illumina Hiseq4000 with approximately 87,000 reads per cell \cite{zheng2017massively}, 4645 tumor-derived single cells from 19 melanoma patient samples, 3114 single cells from four ovarian cancer ascites samples. The second dataset (PBMC-12k merged) contains the 4k PBMCs and the 8k PBMCs\footnote{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k} from a healthy donor, melanoma patient samples, and four ovarian cancer ascites samples. The 8k PBMCs has 8,381 cells detected sequenced on Illumina Hiseq4000 with approximately 92,000 reads per cell \cite{zheng2017massively}.  The performance measure reported by 10-fold cross-validation is the average of the values calculated.

# Load data
SingleRMEA accepts as input a matrix of raw gene counts with genes as rows and cells as columns. The table should have the format shown below.

# Normalization by housekeeping genes 

Importantly, when we merged those different samples, there exist some bath effects. To address this problem, we apply the to select the housekeeping genes for normalization to decrease statistical power. To minimize platform-dependent errors, we select the housekeeping genes for normalization. In this study, we employed the 3804 housekeeping genes (HK_genes.mat) to normalize the single-cell RNA-seq data.

    % restrict to common genes
    [~, ia, ib] = intersect(data1.Properties.RowNames,data2.Properties.RowNames);
    data1 = data1(ia,:);
    data2 = data2(ib,:);
    clear ia ib;

    % load house-keeping genes
    load('HK_genes.mat');

    % find common genes
    [~, ~, ia] = intersect(hk_genes,data1.Properties.RowNames);

    % convert to TPM scale
    data1 = logTrafo(data1,-1);
    data2 = logTrafo(data2,-1);

    % normalize to house-keeping gene expression
    hk_expr = mean([data1{ia,:} data2{ia,:}],1);
    id1 = [true(1,size(data1,2)) false(1,size(data2,2))];
    id2 = [false(1,size(data1,2)) true(1,size(data2,2))];

    data1{:,:} = bsxfun(@times,data1{:,:},mean(hk_expr)./hk_expr(id1));
    data2{:,:} = bsxfun(@times,data2{:,:},mean(hk_expr)./hk_expr(id2));

    % convert back to log scale
    data1 = logTrafo(data1,1);
    data2 = logTrafo(data2,1);

# Dimensionality Reduction (t-SNE)
t-Distributed Stochastic Neighbor Embedding (t-SNE) is a (prize-winning) technique for dimensionality reduction that is particularly well suited for the visualization of high-dimensional datasets. The technique can be implemented via Barnes-Hut approximations, allowing it to be applied on large real-world datasets. 

We obtained the code from the https://lvdmaaten.github.io/tsne/

# Clustering

To design the computational methods for classifying cell types on single-cell RNA-seq datasets, the cell type labels associated with the gene expression profiles are needed to train a machine learning classier. For single-cell RNA-seq dataset, its cell types labels are often unavailable. Therefore, it is necessary to determine the cell type of each cell using the clustering  methods. For clustering single-cell RNA-seq datasets, a multi-step method was developed to divide the cells in single-cell RNA-seq datasets into a number of clusters (also known as cell sub-populations). The multi-step method includes two critical phases: dimensional reduction and clustering. In the first phase, t-Distributed Stochastic Neighbor Embedding (t-SNE) \cite{maaten2008visualizing,van2014accelerating} is applied to identify similar cells. After that, since the cell labels and the number of clusters are available, different clustering algorithms including Density-based spatial clustering of applications with noise (DBSCAN) and CDP are employed for clustering based on the t-SNE mapping on PBMC-4K merged data  (symbol types: squares for PBMC-4K, triangles for the melanoma-data sets, and diamonds for ascites data-sets) as shown in Figure \ref{fig:Figure1}. DBSCAN is
the density associated with a point by counting the number of points in a region of specified radius around the point. The parameters of DBSCAN are set to MinPts=25 and Eps=2 {\color{blue}(the parameters were selected (after some preliminary experiments) so as to result in the clusters generated by the other algorithms used for comparison. However, with different parameters, it is very difficult to ensure that it is the best setting for this single-cell RNA-seq data.)}.

CDP is the clustering algorithm that
efficiently discovers the centers of clusters by finding the density peaks. The parameters of CDP are set to percent = 2. To compare those clustering algorithms, the silhouette coefficient is calculated to measure the performance of cell type label generation by different algorithms. A higher silhouette coefficient score represents a model with well-shaped clusters and vice versa. As depicted in Figure \ref{fig:Figure1}, the CDP can achieve better silhouette coefficient values than DBSCAN. We can see that the t-SNE mapping can be divided into 24 clusters as shown in Figure \ref{fig:Figure2}. Another interesting application of this heatmap is the ability to cluster the cells by their gene expression profiles. Supplementary Figure \ref{fig:FigureS1} summarized the silhouette coefficient scores, the number of clusters on PBMC-12K merged data for DBSCAN and CDP. Those comparison results also show that CDP is very suitable for clustering the PBMC data. Supplementary Figure \ref{fig:FigureS2} shows the heatmap for PBMC-12K merged data.

