RNA-seq data analysis pipeline in R for differential gene expression analysis and functional enrichment, with a focus on visualization, and interpretation of biological data.

DE_and_Enrichment_Analysis.R:

Extracts and pre-processes raw RNA-seq count data and metadata from GEO datasets using GEOquery and DESeq2.
Identifies up- and down-regulated genes using DESeq2, filtering based on adjusted p-values (< 0.01) and log2 fold-change thresholds.
Conducts Gene Ontology (GO) enrichment analysis for biological processes (BP) using clusterProfiler.
Creats plots (PCA, enhanced volcano, heatmap, cnet) for data visualisation

Identifies genes directly regulated by PRC2 (Polycomb Repressive Complex 2) using ChIP-seq overlap data, linking epigenetic modifications (H3K27me3) with RNA-seq expression changes.
