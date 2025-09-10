library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(stringr)

BiocManager::install("vsn")


# TruSeq® Stranded mRNA Library Prep Kit (RS-122-2101, Illumina) was 
# used to prepare these samples. Since this is a stranded library prep
# kit, select from column 3 or 4. In the output table column 4 contains 
# the reads whereas column 3 does not, proceed with using column 4 from 
# each of the samples for downstream analysis.

# Column to extract (4 = Python index 3)
desired_column <- 4  
output_directory <- "STAR_output/"

# Find all ReadsPerGene.out.tab files
files <- list.files(
  output_directory,
  pattern = glob2rx("*.ReadsPerGene.out.tab"),
  full.names = TRUE
)

if(length(files) == 0){
  stop("No ReadsPerGene.out.tab files found in ", output_directory)
}

# Extract sample names from filenames (remove directory and suffix)
sample_names <- str_replace(basename(files), "\\.ReadsPerGene\\.out\\.tab$", "")

# Read files into a list
sample_list <- list()
for(i in seq_along(files)){
  tab <- read.table(files[i], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  if(ncol(tab) < desired_column){
    warning("File ", files[i], " has fewer than ", desired_column, " columns — skipping.")
    next
  }
  
  df <- tab[, c(1, desired_column), drop = FALSE]
  colnames(df) <- c("ensembl_id", sample_names[i])
  sample_list[[ sample_names[i] ]] <- df
  cat("Column", desired_column, "of", sample_names[i], "stored\n")
}

# Check if any valid samples
if(length(sample_list) == 0){
  stop("No valid ReadsPerGene.out.tab files could be processed. Exiting.")
}

# Merge all samples by gene ID
if(length(sample_list) > 1){
  merged <- Reduce(function(x, y) merge(x, y, by = "ensembl_id", all = TRUE), sample_list)
} else {
  merged <- sample_list[[1]]
}

# Split QC and gene counts
qc_rows <- grepl("^N_", merged$ensembl_id)
qc <- merged[qc_rows, , drop = FALSE]
counts <- merged[!qc_rows, , drop = FALSE]

# Write output
if(nrow(counts) > 0){
  write.csv(counts, "raw_counts.csv", row.names = FALSE, quote = FALSE)
  cat("raw_counts.csv written with", nrow(counts), "rows\n")
} else {
  warning("No gene counts to write to raw_counts.csv")
}

if(nrow(qc) > 0){
  write.csv(qc, "qc.csv", row.names = FALSE, quote = FALSE)
  cat("qc.csv written with", nrow(qc), "rows\n")
} else {
  warning("No QC metrics to write to qc.csv")
}

#======================================================================
# Input raw_counts.csv data

data <- read.csv("raw_counts.csv", header = T, row.names = "ensembl_id")
data <- data[,sort(colnames(data))]

head(data)    # check head
colSums(data) # check colum

# ===============================================================================================================
# Differential expression analysis
# Create the DESeq2 object 
# construct the colData, identify biological replicates
condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), 
               rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))

# assign replicates to each sample name to construct colData
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
my_colData


dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)
dds
head(dds@assays@data$counts)  # check the dds object
normalized_counts <- counts(dds, normalized = T) # normalized 
head(normalized_counts)

# Result 
res <- results(dds)
res
res <- results(dds, name = "condition_PC3_Normoxia_vs_LNCAP_Hypoxia")
res <- results(dds, contrast=c("condition","PC3_Normoxia","LNCAP_Hypoxia"))
head(res)

# Log fold change shrinkage for visualization and ranking
library(apeglm)
library(ashr)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_PC3_Normoxia_vs_LNCAP_Hypoxia", type="apeglm")
resLFC

#==================================================================================================================
# Annotation file with GRCh38 with Ensembl ID
annotation <- read.csv("GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
head(annotation)

# Add fourth column
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), 
                                        var = "ensembl_id")
annotated_data <- right_join(annotation, normalized_counts, 
                             by = c("Gene.stable.ID" = "ensembl_id"))
head(annotated_data)

write.csv(annotated_data, file = "gene_annotated_normalized_counts.csv")
vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)
head(assay(vsd), 3)
head(assay(rld), 3)

# ==========================================================================================================
# Visualization
# MA plots
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# condition_PC3_Hypoxia_vs_LNCAP_Hypoxia
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_PC3_Hypoxia_vs_LNCAP_Hypoxia", type="apeglm")
resLFC
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

## because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

png("MAplots_bulk.png", width=3000, height=1200, res=300)
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

#=============================================================================================================
# Distance Plot
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors) 
}
plotDists(vsd)

#--------------------------------------------------------------------
# Variable Genes Heatmap
variable_gene_heatmap <- function(vsd, num_genes = 500, annotation, title = "") {
  # Colors
  colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)[256:1]
  
  # Top variable genes
  mat <- assay(vsd)
  top_genes <- mat[order(rowVars(mat), decreasing = TRUE)[1:num_genes], ]
  top_genes <- top_genes - rowMeans(top_genes, na.rm = TRUE)
  
  # Replace IDs with gene names
  rownames(top_genes) <- annotation$Gene.name[
    match(rownames(top_genes), annotation$Gene.stable.ID)
  ]
  
  # Sample annotations
  ann <- as.data.frame(colData(vsd))
  ann$sizeFactor <- NULL
  
  # Heatmap
  pheatmap::pheatmap(top_genes, color = colors, annotation_col = ann,
                     fontsize_col = 8, fontsize_row = 250/num_genes,
                     border_color = NA, main = title)
}

variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation, title = "Variable Genes")

#----------------------------------------------------------------------
# PCA Plot
png("PCA_New.png", width=2000, height=1200, res=300)
plot_PCA <- function(vsd.obj) {
  d <- plotPCA(vsd.obj, intgroup="condition", returnData=TRUE)
  v <- round(100 * attr(d, "percentVar"))
  ggplot(d, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    ggrepel::geom_text_repel(aes(label=name), color="black", size=3) +
    labs(x=paste0("PC1 (",v[1],"%)"), y=paste0("PC2 (",v[2],"%)"))
}
print(plot_PCA(vsd))
dev.off()

plot_PCA(vsd)

#-----------------------------------------------------------
# PlotCounts
library(degPlot)

plotCounts(dds, gene="ENSG00000146678", intgroup="condition")
d <- plotCounts(dds, gene="ENSG00000146678", intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0), aes (color = condition)) + 
  scale_y_log10(breaks=c(25,100,400))

head(annotation)
head(dds)

#-Annotated plot counts
# Load annotation table
annotation <- read.csv("GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
head(annotation)


# Select gene of interest by gene symbol
gene_name <- "IGFBP1"
gene_id   <- annotation$Gene.stable.ID[annotation$Gene.name == gene_name]

# Get normalized counts for that gene
d <- plotCounts(dds, gene=gene_id, intgroup="condition", returnData=TRUE)

# Plot
ggplot(d, aes(x=condition, y=count, color=condition)) + 
  geom_point(position=position_jitter(w=0.1, h=0), size=3) + 
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,5)) +
  labs(
    title = paste0("Expression of ", gene_name, " (", gene_id, ")"),
    x = "Condition",
    y = "Normalized count"
  ) +
  theme_classic(base_size=14) +   # cleaner look
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1)  # add box
  )

#==================================================================================================================
# Making subsets -------LNCAP and PC3--------------------------
generate_DESeq_object <- function(my_data, groups) {
  # Select samples for the two groups
  selected_cols <- grep(paste0("^(", paste(groups, collapse = "|"), ")"), colnames(my_data))
  count_data <- my_data[, selected_cols]
  
  # Create condition labels
  condition <- factor(ifelse(grepl(paste0("^", groups[1]), colnames(count_data)),
                             groups[1], groups[2]))
  
  col_data <- data.frame(condition = condition,
                         row.names = colnames(count_data))
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = col_data,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  return(dds)
}

LNCAP <- generate_DESeq_object(data, c("LNCAP_Hypoxia", "LNCAP_Normoxia"))
PC3 <- generate_DESeq_object(data, c("PC3_Hypoxia", "PC3_Normoxia"))
head(LNCAP)

#===================================================================================================================
#Extracting DE results
results(LNCAP, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))
results(dds, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))

library(DESeq2)
library(dplyr)

# Get DESeq2 results
res <- results(dds, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)
lncap_counts <- norm_counts[, grepl("LNCAP", colnames(norm_counts))]

# Calculate average CPM (counts per million)
cpm <- edgeR::cpm(norm_counts)  # from edgeR
avg_cpm <- rowMeans(cpm[, grepl("LNCAP", colnames(cpm))])
#*******************************************************************************

generate_DE_results <- function(dds, groups, annotation) {

# DESeq2 results
  res <- results(dds, contrast = c("condition", groups[1], groups[2]))
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)
  
# Remove unwanted DESeq2 columns
  res_df <- res_df[, !(colnames(res_df) %in% c("stat", "lfcSE"))]
  
# Normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  cpm <- edgeR::cpm(norm_counts)  
  avg_cpm <- rowMeans(cpm[, grepl(groups[1], colnames(cpm)) | grepl(groups[2], colnames(cpm))])
  
  counts_df <- as.data.frame(norm_counts[, grepl(groups[1], colnames(norm_counts)) | grepl(groups[2], colnames(norm_counts))])
  counts_df$avg_cpm <- avg_cpm
  counts_df$ensembl_id <- rownames(counts_df)
  
# Join everything
  combined <- res_df %>%
    left_join(annotation, by = c("ensembl_id" = "Gene.stable.ID")) %>%
    left_join(counts_df, by = "ensembl_id") %>%
    select(ensembl_id, baseMean, log2FoldChange, pvalue, padj, 
           Gene.name, Gene.type, avg_cpm, everything())
  
# Summary
  total_genes <- nrow(combined)
  sig_genes   <- sum(combined$padj < 0.001, na.rm = TRUE)
  lfc_genes   <- sum(combined$padj < 0.001 & abs(combined$log2FoldChange) > 0.5, na.rm = TRUE)
  cpm_genes   <- sum(combined$padj < 0.001 & combined$avg_cpm > 2, na.rm = TRUE)
  
  message(sprintf("## For %s_vs_%s: %d total genes", groups[1], groups[2], total_genes))
  message(sprintf("## %d padj < 0.001", sig_genes))
  message(sprintf("## %d padj < 0.001 & |log2FC| > 0.5", lfc_genes))
  message(sprintf("## %d padj < 0.001 & avg CPM > 2", cpm_genes))
  
  return(combined)
}

# Run and save CSV
lncap_output <- generate_DE_results(dds, c("LNCAP_Hypoxia", "LNCAP_Normoxia"), annotation)
write.csv(lncap_output, file = "LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv", row.names = FALSE)

# Read back
res_LNCAP <- read.csv("LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv", header = TRUE)
head(res_LNCAP)

pc3_output <- generate_DE_results(dds, c("PC3_Hypoxia", "PC3_Normoxia"), annotation)
write.csv(pc3_output, file = "PC3_Hypoxia_vs_PC3_Normoxia_allgenes.csv",
          row.names = FALSE)

res_PC3 <- read.csv("PC3_Hypoxia_vs_PC3_Normoxia_allgenes.csv", header = TRUE)
head(res_PC3)

#==============================================================================================================
# Differential gene heatmap
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(grid)

# Load CSV and matrix
res_LNCAP <- read.csv("LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv",
                stringsAsFactors = FALSE, check.names = FALSE)
norm_counts <- as.matrix(res_LNCAP[, c("LNCAP_Hypoxia_S1","LNCAP_Hypoxia_S2",
                                 "LNCAP_Normoxia_S1","LNCAP_Normoxia_S2")])
rownames(norm_counts) <- res_LNCAP$ensembl_id

# Heatmap of top upregulated genes
sig <- res_LNCAP %>%
  filter(!is.na(Gene.name), Gene.name != "",
         !is.na(log2FoldChange), padj < 1e-3, log2FoldChange > 0) %>%
  slice_max(log2FoldChange, n = 30)

mat <- t(scale(t(norm_counts[sig$ensembl_id, , drop = FALSE])))
mat[is.na(mat)] <- 0
rownames(mat) <- sig$Gene.name

png("LNCAP_upregulated_heatmap.png", width = 800, height = 600)
grid::grid.draw(pheatmap(mat,
                         color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(256)),
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         fontsize_col = 16, fontsize_row = 12,
                         border_color = NA)
                $gtable)
dev.off()

#--------------PC3-------------------------------------------
# Load CSV and matrix
res_PC3 <- read.csv("PC3_Hypoxia_vs_PC3_Normoxia_allgenes.csv",
                      stringsAsFactors = FALSE, check.names = FALSE)
norm_counts <- as.matrix(res_PC3[, c("PC3_Hypoxia_S1","PC3_Hypoxia_S2",
                                       "PC3_Normoxia_S1","PC3_Normoxia_S2")])
rownames(norm_counts) <- res_PC3$ensembl_id

# Heatmap of top upregulated genes
sig <- res_PC3 %>%
  filter(!is.na(Gene.name), Gene.name != "",
         !is.na(log2FoldChange), padj < 1e-3, log2FoldChange > 0) %>%
  slice_max(log2FoldChange, n = 30)

mat <- t(scale(t(norm_counts[sig$ensembl_id, , drop = FALSE])))
mat[is.na(mat)] <- 0
rownames(mat) <- sig$Gene.name

png("PC3_upregulated_heatmap.png", width = 800, height = 600)
grid::grid.draw(pheatmap(mat,
                         color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(256)),
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         fontsize_col = 16, fontsize_row = 12,
                         border_color = NA)
                $gtable)
dev.off()
#=====================================================================================================================
# Volcano plot-
library(dplyr)
library(ggplot2)
library(ggrepel)

res_LNCAP <- read.csv("LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv")

plot_volcano <- function(res_LNCAP, padj_cutoff, nlabel = 10, label.by = c("padj", "log2FoldChange")) {
  label.by <- match.arg(label.by)
  
  res_LNCAP <- res_LNCAP %>%
    mutate(significance = ifelse(padj < padj_cutoff, "significant", "not significant")) %>%
    filter(!is.na(significance))
  
  sig <- filter(res_LNCAP, significance == "significant")
  
  # Select top and bottom genes based on label.by
  if(label.by == "padj") {
    top <- sig %>% arrange(padj) %>% slice_head(n = nlabel)
    bottom <- sig %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% slice_head(n = nlabel)
  } else {
    top <- sig %>% arrange(desc(log2FoldChange)) %>% slice_head(n = nlabel)
    bottom <- sig %>% arrange(log2FoldChange) %>% slice_head(n = nlabel)
  }
  
  ggplot(res_LNCAP, aes(log2FoldChange, -log10(padj), color = significance)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("significant" = "red", "not significant" = "black")) +
    geom_text_repel(data = top, aes(label = Gene.name), size = 6, color = "red", max.overlaps = 50) +
    geom_text_repel(data = bottom, aes(label = Gene.name), size = 6, color = "#619CFF", max.overlaps = 50) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(x = "Log2 Fold Change", y = "-log10(padj)") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_text(face = "bold", size = 25),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.position = c(0.94, 0.94)
    )
}
# Create the plot view 
volcano_plot <- plot_volcano(res, 0.0005, nlabel = 20, label.by = "padj")
volcano_plot

# Create and save the plot
ggsave("volcano_plot.png", 
       plot = plot_volcano(res, 0.0005, nlabel = 20, label.by = "padj"),
       width = 18,       # Increase width
       height = 10,      # Increase height
       dpi = 300,        # Higher resolution
       bg = "white")

#---------PC3----------------------------------------------------------------
# Volcano plot-
library(dplyr)
library(ggplot2)
library(ggrepel)

res_PC3 <- read.csv("PC3_Hypoxia_vs_PC3_Normoxia_allgenes.csv")

plot_volcano <- function(res_PC3, padj_cutoff, nlabel = 10, label.by = c("padj", "log2FoldChange")) {
  label.by <- match.arg(label.by)
  
  res_PC3 <- res_PC3 %>%
    mutate(significance = ifelse(padj < padj_cutoff, "significant", "not significant")) %>%
    filter(!is.na(significance))
  
  sig <- filter(res_PC3, significance == "significant")
  
  # Select top and bottom genes based on label.by
  if(label.by == "padj") {
    top <- sig %>% arrange(padj) %>% slice_head(n = nlabel)
    bottom <- sig %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% slice_head(n = nlabel)
  } else {
    top <- sig %>% arrange(desc(log2FoldChange)) %>% slice_head(n = nlabel)
    bottom <- sig %>% arrange(log2FoldChange) %>% slice_head(n = nlabel)
  }
  
  ggplot(res_PC3, aes(log2FoldChange, -log10(padj), color = significance)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("significant" = "red", "not significant" = "black")) +
    geom_text_repel(data = top, aes(label = Gene.name), size = 6, color = "red", max.overlaps = 50) +
    geom_text_repel(data = bottom, aes(label = Gene.name), size = 6, color = "#619CFF", max.overlaps = 50) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(x = "Log2 Fold Change", y = "-log10(padj)") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_text(face = "bold", size = 25),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.position = c(0.94, 0.94)
    )
}
# Create the plot view 
volcano_plot <- plot_volcano(res_PC3, 0.0005, nlabel = 20, label.by = "padj")
volcano_plot

# Create and save the plot
ggsave("volcano_plot_PC3.png", 
       plot = plot_volcano(res_PC3, 0.0005, nlabel = 20, label.by = "padj"),
       width = 18,       # Increase width
       height = 12,      # Increase height
       dpi = 300,        # Higher resolution
       bg = "white")

# ===============================================================================================================
# Making subsets -------LNCAP and PC3--------------------------
generate_DESeq_object <- function(my_data, groups) {

  # Select samples for the two groups
  selected_cols <- grep(paste0("^(", paste(groups, collapse = "|"), ")"), colnames(my_data))
  count_data <- my_data[, selected_cols]
  
# Create condition labels
  condition <- factor(ifelse(grepl(paste0("^", groups[1]), colnames(count_data)),
                             groups[1], groups[2]))
  
  col_data <- data.frame(condition = condition,
                         row.names = colnames(count_data))
  
# Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = col_data,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  return(dds)
}

LNCAP <- generate_DESeq_object(data, c("LNCAP_Hypoxia", "LNCAP_Normoxia"))
PC3 <- generate_DESeq_object(data, c("PC3_Hypoxia", "PC3_Normoxia"))
head(LNCAP)

library(gridExtra)

lncap_vsd <- vst(lncap, blind = T)
pc3_vsd <- vst(pc3, blind = T)
LNCAP_vis <- variable_gene_heatmap(lncap_vsd, 30, annotation = annotation, 
                           title = "LNCaP variable genes")
PC3_vis <- variable_gene_heatmap(pc3_vsd, 30, annotation = annotation, 
                           title = "PC3 variable genes")

grid.arrange(LNCAP_vis$gtable, PC3_vis$gtable, ncol = 2)

# ==================================================================================================================
# Gene set enrichment analysis Hall mark_pathway.gmt

library(fgsea)
# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("h.all.v7.0.symbols.gmt")
head(names(hallmark_pathway))

# Extract ranked genes------------
library(dplyr)
library(stringr)

res_prot <- res_LNCAP %>%
  filter(Gene.type == "protein_coding") %>%
  select(Gene.name, log2FoldChange) %>%
  filter(!is.na(Gene.name), Gene.name != "") %>%
  na.omit() %>%
  mutate(Gene.name = str_to_upper(str_replace_all(Gene.name, "[[:space:]]+", "_")))

# Write sorted .rnk with original header
write.table(
  res_prot[order(-res_prot$log2FoldChange), ],
  file = "LNCAP_Hypoxia_vs_LNCAP_Normoxia_proteins.rnk",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# load the ranked list
lncap_ranked_list <- read.table("LNCAP_Hypoxia_vs_LNCAP_Normoxia_proteins.rnk", header = T, stringsAsFactors = F)
head(lncap_ranked_list)

# Ranked list
prepare_ranked_list <- function(df) {
  df %>%
    group_by(Gene.name) %>%
    summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(log2FoldChange)) %>%
    drop_na() %>%
    deframe()
}

lncap_ranked_list <- prepare_ranked_list(lncap_ranked_list)
head(lncap_ranked_list)

# halmark_pathways
fgsea_results_PW <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = lncap_ranked_list,
                                 minSize = 15,
                                 maxSize = 500)

fgsea_results_PW %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

# Waterfall Plot
library(dplyr)
library(ggplot2)
library(stringr)

waterfall_plot <- function(fgsea_results_PW, title = "Waterfall Plot") {
  fgsea_results_PW %>%
    mutate(pathway = str_remove(pathway, "^HALLMARK_")) %>%
    ggplot(aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = title) +
    theme(axis.text.y = element_text(size = 30),       # increase y-axis labels
          axis.text.x = element_text(size = 30),       # increase x-axis labels
          axis.title = element_text(size = 45, face = "bold"),  # axis titles
          plot.title = element_text(hjust = 0.5, size = 60, face = "bold"), # title
          legend.text = element_text(size = 20),      # legend text
          legend.title = element_text(size = 45, face = "bold"))
}

# Display plot
# waterfall_plot(fgsea_results_PW, "Hallmark pathways altered by hypoxia in LNCaP cells")

# Save the plot
ggsave("LNCAP_waterfall_PW.png", waterfall_plot(fgsea_results_PW, "Hallmark pathways altered by hypoxia in LNCaP cells"),
       width = 30, height = max(6, nrow(fgsea_results_PW)/2))

#===============================================================================
# Using msigdbr for analysis
library(msigdbr)
library(dplyr)
library(tibble)

# Get all hallmark gene sets for human
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
# Filter for hypoxia only
hallmark_hypoxia_genes <- subset(hallmark, gs_name == "HALLMARK_HYPOXIA")
# Extract just the gene symbols
hypoxia_genes <- unique(hallmark_hypoxia_genes$gene_symbol)
head(hypoxia_genes)

#-------------------------------
# LNCAP ranked genes
lncap_ranked_list <- read.table("LNCAP_Hypoxia_vs_LNCAP_Normoxia_proteins.rnk", header = T, stringsAsFactors = F)
head(lncap_ranked_list)

# hallmark_pathways list exactly like gmtPathways() 
hallmark.pathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()
head(names(hallmark.pathways))

# format ranked list for the fgsea() function

prepare_ranked_list <- function(df) {
  df %>%
    group_by(Gene.name) %>%
    summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(log2FoldChange)) %>%
    drop_na() %>%
    deframe()
}

lncap_ranked_list <- prepare_ranked_list(lncap_ranked_list)
head(lncap_ranked_list)

# halmark_pathways
fgsea_results_msigdbr <- fgseaMultilevel(pathways = hallmark.pathways,
                                    stats = lncap_ranked_list,
                                    minSize = 15,
                                    maxSize = 500)

fgsea_results_msigdbr %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

# Waterfall Plot
library(dplyr)
library(ggplot2)
library(stringr)

waterfall_plot <- function(fgsea_results_msigdbr, title = "Waterfall Plot") {
  fgsea_results_msigdbr %>%
    mutate(pathway = str_remove(pathway, "^HALLMARK_")) %>%
    ggplot(aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = title) +
    theme(axis.text.y = element_text(size = 30),       # increase y-axis labels
          axis.text.x = element_text(size = 30),       # increase x-axis labels
          axis.title = element_text(size = 45, face = "bold"),  # axis titles
          plot.title = element_text(hjust = 0.5, size = 60, face = "bold"), # title
          legend.text = element_text(size = 20),      # legend text
          legend.title = element_text(size = 45, face = "bold"))
}

# Display plot
# waterfall_plot(fgsea_results_PW, "Hallmark pathways altered by hypoxia in LNCaP cells")

# Save the plot
ggsave("LNCAP_waterfall_msigdbr.png", waterfall_plot(fgsea_results_PW, "Hallmark pathways altered by hypoxia in LNCaP cells"),
       width = 30, height = max(6, nrow(fgsea_results_PW)/2))


# Suppose your DE results are in a data frame 'res_prot' with columns Gene.name and log2FoldChange
lncap_ranked_list <- res_prot$log2FoldChange        # numeric vector
names(lncap_ranked_list) <- res_prot$Gene.name     # assign gene names
lncap_ranked_list <- sort(lncap_ranked_list, decreasing = TRUE)

hypoxia_genes <- list(Hypoxia = hypoxia_genes)

#=====================================================================================================================
# Enriched specific pathways 
plot_enrichment <- function(geneset, pathway, ranked_list) {
  enrich_plot <- plotEnrichment(geneset[[pathway]], ranked_list)
  enrich_plot$layers[[1]]$aes_params$colour <- "blue"
  enrich_plot + labs(title = pathway)
}

# Example usage
plot_enrichment(hallmark_pathway, "HALLMARK_HYPOXIA", lncap_ranked_list)

ggsave("HALLMARK_HYPOXIA_enriched.png",
       plot_enrichment(hallmark_pathway, "HALLMARK_HYPOXIA", lncap_ranked_list),
       width = 8,
       height = 4)


# Negatively enriched pathway 
plot_enrichment(hallmark_pathway, "HALLMARK_PEROXISOME",
                lncap_ranked_list)

ggsave("HALLMARK_PEROXISOME_negatively_enriched.png",
       plot_enrichment(hallmark_pathway, "HALLMARK_PEROXISOME",
                       lncap_ranked_list),
       width = 8,
       height = 4)
