################
#### PART 0 ####
################
#### --------------------- INSTALLING AND LOADING PACKAGES -------------------------------

library(DESeq2)
library(tidyverse)
library(GEOquery)
library(here)
library(pheatmap)
library(biomaRt)
#BiocManager::install("biomaRt")
library(EnhancedVolcano)
#BiocManager::install("EnhancedVolcano")
library(org.Mm.eg.db) # Annotation library for mouse
#BiocManager::install("org.Mm.eg.db")
library(clusterProfiler) # Enrichment package
#BiocManager::install("clusterProfiler")
library(enrichplot)

################
#### PART 1 ####
################
#### --------------------- GETTING COLDATA, ASSAYDATA, ROWDATA FOR DESEQ ------------------------------- 

counts <- read_tsv(here("data" , "GSE169312_ReadCount.txt.gz"))

### GETTING ASSAY DATA
assayData <- as.data.frame(counts[,c(2:7)])
rownames(assayData) <- counts$ID
assayData <- assayData |> 
  relocate(colnames(assayData)[c(4:6)])
# in the assay data, EED knock-out samples come before the wild types, but in the metadata, the wild types come first
# bring forward the wildtype columns to the front

# Convert all the counts in AsssayData to integers
assayData <- as.data.frame(sapply(assayData, as.integer))
rownames(assayData) <- counts$ID
names(assayData) <- c("WT1", "WT2", "WT3", "cKO1", "cKO2", "cKO3")

### GETTING ROW DATA
rowData <- as.data.frame(counts[,c(1)])
rownames(rowData) <- counts$ID

### GETTING COL DATA / DESCRIBING THE EXPERIMENT 
gse169312 <- getGEO("GSE169312")[[1]]
colData <- pData(gse169312) |> janitor::clean_names()
names(colData) <- names(colData) |> 
  stringr::str_remove("_ch1")
# Cleans the names and removes "_ch1"

# Since we are only interested in RNA-seq data for DESeq, exclude the information regarding the ChIP-seq samples.
colData <- colData[1:6,]

# Changing the genotype names, make the "genotype" variable a factor.
colData$genotype <- as.factor(colData$genotype)
colData$genotype <- factor(rep(c("WT", "cKO"), each = 3))
colData$genotype <- forcats::fct_relevel(colData$genotype, "WT")
colData$genotype

# Rename the row names of col data
rownames(colData) <- c("WT1", "WT2", "WT3", "cKO1", "cKO2", "cKO3")

# check if the col names of assayData is equal to rownames of colData (metadata). 
all(rownames(colData)==colnames(assayData))

################
#### PART 2 ####
################
#### ---------------------- PEFORMING DESEQ --------------------- ####

# Getting summarized data
gse169312.se <- SummarizedExperiment(assays = assayData,
                                     colData = colData,
                                     rowData = rowData)
#gse169312.se

# Building DEseq Data Set
gse169312.dds <- DESeqDataSetFromMatrix(countData = assayData,
                                        colData = colData,
                                        rowData = rowData,
                                        design = ~ genotype)

# Running DE seq on count data, and extracting results
gse169312.dds <- DESeq(gse169312.dds)
res <- results(gse169312.dds)

################
#### PART 3 ####
################
#### ---------------------- SEPERATING UP- AND DOWN- REGULATED GENES --------------------- ####

# View the number of reads that have a non-zero total read count, and other summary statistics 
summary(res)
# Quick visualization of how the results would look like; results are not filtered using any arbitrary threshold yet, nor are they sorted.
res

# Finding significant genes (padj value <0.01) 
# remove reads with NA padj values (which were reponsible for ~60% of the records)
res.sig <- res[(res$padj<0.01)&(!is.na(res$padj)),]
#res.sig

# Finding UP regulated genes (LFC >= 1) and ordering them by padj.
res.sig.up <- res.sig[res.sig$log2FoldChange >= 1,]
res.sig.up <- res.sig.up[order(res.sig.up$padj,
                               decreasing=F),]
#res.sig.up
# Preliminary results shows that arbitrary thresholds of padj < 0.01 and LFC > 1 yielded a total of 649 up-regulated genes in the EED-knockout samples.

# Saving the Ensembl gene IDs (without version) of significantly up-regulated genes into a .txt file
res.sig.up.export <- data.frame(rownames(res.sig.up))
write.table(sapply(res.sig.up.export$rownames.res.sig.up.,
                   function(x) strsplit(x,"[.]")[[1]][1]), 
            here("results" , "Group10_sig_up_genes.txt"), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)

# Finding DOWN regulated genes (LFC <= 1) and ordering them by padj.
res.sig.down <- res.sig[res.sig$log2FoldChange <= (-1),]
res.sig.down <- res.sig.down[order(res.sig.down$padj,
                                   decreasing=F),]
#res.sig.down
# Preliminary results shows that arbitrary thresholds of padj < 0.01 and LFC < 1 yielded a total of 775 down-regulated genes in the EED-knockout samples.
# Saving the Ensembl gene IDs (without version) of significantly down-regulated genes into a .txt file
res.sig.down.export <- data.frame(rownames(res.sig.down))
write.table(sapply(res.sig.down.export$rownames.res.sig.down.,
                   function(x) strsplit(x,"[.]")[[1]][1]), 
            here("results" , "Group10_sig_down_genes.txt"), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)

################
#### PART 4 ####
################
#### -------------------- OBTAINING PLOTS ------------------------------- ####

# PLOT 1: PCA plot
plotPCA(rlog(gse169312.dds), intgroup="genotype") +
  theme_minimal() + 
  labs(x = "PC1: 74% explained variance",
       y = "PC2: 9% explained variance")

# PLOT 2: Enhanced volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = c(rownames(res.sig.up)[1:10],rownames(res.sig.down)[1:10]),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 2.5) + 
  labs(x = "log2FoldChange",
       y = "-log10(padj)",
       title = "Volcano Plot for Differentially Expressed Genes",
       subtitle = "Total of 55,385 genes",
       caption = "",
       color = "Significance") +
  theme_minimal() +
  xlim(-12,12) + 
  theme(legend.position = "top")

# clustered heat map
rlog.norm.counts <- assay(rlog(gse169312.dds))

# PLOT 3: HEATMAP FOR ALL UP- AND DOWN-REGULATED GENES WITHOUT ENSEMBL GENE IDs (too many genes)
norm.mat_all <- rlog.norm.counts[c(rownames(res.sig.up),
                                   rownames(res.sig.down)),]
pheatmap(norm.mat_all,scale="row",show_rownames = FALSE)

# PLOT 4: HEATMAP FOR TOP 25 UP- AND DOWN-REGULATED GENES WITH ENSEMBL GENE IDs
norm.mat <- rlog.norm.counts[c(rownames(res.sig.up)[1:25],
                               rownames(res.sig.down)[1:25]),]
#pheatmap(norm.mat,scale="row")

# Replacing Emsembl gene symbols with official gene symbols for clustered heatmap of top/bottom 25 genes
diffgenes <- data.frame(ensembl_gene_id_version=rownames(norm.mat),
                        ensembl_gene_id=NA)
diffgenes$ensembl_gene_id <- sapply(diffgenes$ensembl_gene_id_version,
                                    function(x)
                                      strsplit(x,"[.]")[[1]][1])
diffgenes
# getting rid of the version

# (two possible commands) recieve the mouse gene annotation data from Ensembl
#mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
mart <- useEnsembl("ensembl", "mmusculus_gene_ensembl")

z <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
           filters = "ensembl_gene_id",
           values = diffgenes$ensembl_gene_id,
           mart = mart)

diffgenes.symbol <- merge(diffgenes, z)

symbols <- diffgenes.symbol[match(rownames(norm.mat),
                                  diffgenes.symbol$ensembl_gene_id_version),"mgi_symbol"]

rownames(norm.mat) <- symbols

pheatmap(norm.mat,scale="row",)


################
#### PART 5 ####
################
#### ----------------------  ENRICHMENT ANALYSIS OF UP/DOWN REGULATED GENES ------------------------------ ####


# DOING ENRICHMENT ANALYSIS FOR ALL UPREGULATED GENES
DE_genes_up <- c(rownames(res.sig.up))

# version number from DE_genes_up
DE_genes_up_noversion <- data.frame(ensembl_gene_id_version=DE_genes_up,
                                    ensembl_gene_id=NA)
DE_genes_up_noversion$ensembl_gene_id <- sapply(DE_genes_up_noversion$ensembl_gene_id_version,
                                                function(x)
                                                  strsplit(x,"[.]")[[1]][1])
#DE_genes_up_noversion$ensembl_gene_id

# Need to remove version for universe genes (universe will be the ensmble gene IDS of ALL genes)

universe_noversion <- data.frame(ensembl_gene_id_version=rownames(res),
                                 ensembl_gene_id=NA)
universe_noversion$ensembl_gene_id <- sapply(universe_noversion$ensembl_gene_id_version,
                                             function(x)
                                               strsplit(x,"[.]")[[1]][1])
#universe_noversion$ensembl_gene_id

# Gene ontology enrichment analysis for up-regulated genes
ego_up <- enrichGO(gene          = DE_genes_up_noversion$ensembl_gene_id,
                   keyType       = "ENSEMBL",
                   universe      = universe_noversion$ensembl_gene_id,
                   OrgDb         = org.Mm.eg.db,
                   # Fischer's exact test would be done on all categories of BP/CC/MF.
                   ont           = "BP",  # Can be BP, CC, and MF.
                   # Due to the large number of Fischer's exact tests performed
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   # FDR (p-adjusted) cut off
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

View(as_tibble(ego_up))

# Enrich plot
barplot(ego_up,showCategory=12)

# Dotplot
dotplot(ego_up) #+ theme_minimal() 

#cnetplot
cnetplot(ego_up,showCategory=5,cex_label_gene=0.5)

# emap plot (semantic similarity, whereby the different categories are grouped graphically according to their properties)
emapplot(pairwise_termsim(ego_up),showCategory=30,cex_category=0.5)

# DOING ENRICHMENT ANALYSIS FOR DOWN UPREGULATED GENES


DE_genes_down <- c(rownames(res.sig.down))

# need to remove version number from DE_genes_up
DE_genes_down_noversion <- data.frame(ensembl_gene_id_version=DE_genes_down,
                                      ensembl_gene_id=NA)
DE_genes_down_noversion$ensembl_gene_id <- sapply(DE_genes_down_noversion$ensembl_gene_id_version,
                                                  function(x)
                                                    strsplit(x,"[.]")[[1]][1])
#DE_genes_up_noversion$ensembl_gene_id

# universe genes will be the same as before
# Gene ontology enrichment analysis for down-regulated genes
ego_down <- enrichGO(gene          = DE_genes_down_noversion$ensembl_gene_id,
                     keyType       = "ENSEMBL",
                     universe      = universe_noversion$ensembl_gene_id,
                     OrgDb         = org.Mm.eg.db,
                     # Fischer's exact test would be done on all categories of BP/CC/MF.
                     ont           = "BP",  # Can be BP, CC, and MF.
                     # Due to the large number of Fischer's exact tests performed
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     # FDR (p-adjusted) cut off
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

View(as_tibble(ego_down))

# Enrich plot
barplot(ego_down,showCategory=12)

# Dotplot
dotplot(ego_down)

#cnetplot
cnetplot(ego_down,showCategory=5,cex_label_gene=0.5)

# emap plot
# semantic similarity, whereby the different categories are grouped graphically according to 
# their properties
emapplot(pairwise_termsim(ego_down),showCategory=25,cex_category=0.5)





################
#### PART 6 ####
################
#### ---------------------- OBTAINING GENES DIRECTLY REGULATED BY PRC2 ------------------------------ ####

# Read in differentially expressed genes that were upregulated - saved as a .txt file from before
up_genes <- read_table(here("results", "Group10_sig_up_genes.txt"), col_names = FALSE)
#up_genes <- read_table("../results/Group10_sig_up_genes.txt", col_names = FALSE)

# Read in the provided tsv file which contained information of 
# overlaps of ChIP-seq peaks with every single gene in the mouse genome (T/F)
genes_regulated_by_PRC2 <- read_tsv(here("results", "gene_overlaps.tsv"))|>
  # Filter for only up-regulated genes
  dplyr::filter(gene_id %in% up_genes$X1) |>
  # Filter for genes that had H3K27me3 peaks present 
  # for wild type, but not for EED cKO
  # These are genes that are directly regulated by EED/PRC2 complex, as the loss of PRC2
  # function in EED cKO would result in loss of H3K27me3 markers in EED cKO mice as well
  dplyr::filter(WT_H3K27me3 == TRUE, KO_H3K27me3 == FALSE) |>
  # Filter for genes which acetylation markers remained the same regardless of genotype
  # since we are only looking for genes "directly regulated" by trimethylase PRC2.
  dplyr::filter((WT_H3K27ac == FALSE & KO_H3K27ac == FALSE) | (WT_H3K27ac == TRUE & KO_H3K27ac == TRUE))

# Save into a .txt file for keepsake and later use in PART 7 - 
write.table(genes_regulated_by_PRC2$gene_id, here("results", "genes_regulated_by_PRC22.txt"), 
            quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

#write.table(genes_regulated_by_PRC2$gene_id, "results/genes_regulated_by_PRC22.txt", 
#            quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)


################
#### PART 7 ####
################
#### ----------------------  ENRICHMENT ANALYSIS OF GENES DIRECTLY REGULATED BY PRC2 ------------------------------ ####

genes_regulated_by_PRC2$gene_id

# no need to remove version number here
# universe genes will be the same as before

# Gene ontology enrichment analysis for genes directly regulated by PRC2

ego_direct <- enrichGO(gene          = genes_regulated_by_PRC2$gene_id,
                       keyType       = "ENSEMBL",
                       universe      = universe_noversion$ensembl_gene_id,
                       OrgDb         = org.Mm.eg.db,
                       # Fischer's exact test would be done on all categories of BP/CC/MF.
                       ont           = "BP",  # Can be BP, CC, and MF.
                       # Due to the large number of Fischer's exact tests performed
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       # FDR (p-adjusted) cut off
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

View(as_tibble(ego_direct))

# Enrich plot
barplot(ego_direct,showCategory=12)

# Dotplot
dotplot(ego_direct)

#cnetplot
cnetplot(ego_direct,showCategory=5,cex_label_gene=0.5)

# emap plot
# semantic similarity, whereby the different categories are grouped graphically according to 
# their properties
emapplot(pairwise_termsim(ego_direct),showCategory=25,cex_category=0.5)


# ------ ALL ENRICHMENT ANALYSIS PLOTS ( for easy viewing) -----

# ego_up = all up regulated genes
# ego_down = all down regulated genes
# ego_direct = genes directly affected by the EED ko

barplot(ego_up,showCategory=12)
barplot(ego_down,showCategory=12)
barplot(ego_direct,showCategory=12)

dotplot(ego_up)
dotplot(ego_down)
dotplot(ego_direct)

cnetplot(ego_up,showCategory=5,cex_label_gene=0.5)
cnetplot(ego_down,showCategory=5,cex_label_gene=0.5)
cnetplot(ego_direct,showCategory=5,cex_label_gene=0.5)

emapplot(pairwise_termsim(ego_up),showCategory=25,cex_category=0.5)
emapplot(pairwise_termsim(ego_down),showCategory=25,cex_category=0.5)
emapplot(pairwise_termsim(ego_direct),showCategory=25,cex_category=0.5)
