# Set working directory
setwd('~/Documents/QMUL/Crohns')

# Load libraries, register cores
library(data.table)
library(tximport)
library(DESeq2)
library(IHW)
library(tidyverse)
library(bioplotr)
library(BiocParallel)
register(MulticoreParam(8))

# Import data
clin <- fread('./Data/clinical.csv')
t2g <- readRDS('./Data/Hs91.t2g.rds')
anno <- fread('./Data/Hs.anno.csv')
rld <- readRDS('./Data/rld.rds')

# Remove outliers
outliers <- c('1911B-D', '2502A-F', '2502A-CH')
clin <- clin[!SampleID %in% outliers]

# Run tximport
files <- file.path('./Data/Kallisto', clin$SampleID, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)

# Build DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData = clin, 
                                design = ~ Patient + Condition)

# Run DESeq2
dds <- DESeq(dds, parallel = TRUE)

# Prepare gene index and gene name table
idx <- rownames(dds)
e2g <- t2g %>% 
  select(gene_id, gene_name) %>%
  distinct(.)

# Output function
output <- function(trt) {
  lfc <- lfcShrink(dds, contrast = c('Condition', trt, 'Control'),
                   parallel = TRUE, quiet = TRUE) %>%
    as_tibble(.) %>%
    mutate(gene_id = idx) %>%
    rename(logFC = log2FoldChange) %>%
    select(gene_id, logFC)
  res <- results(dds, contrast = c('Condition', trt, 'Control'), tidy = TRUE, 
                 filterFun = ihw, alpha = 0.05) %>%
    na.omit(.) %>%
    mutate(AvgExpr = log2(baseMean)) %>%
    rename(gene_id = row, 
           p.value = pvalue,
           q.value = padj) %>%
    inner_join(lfc, by = 'gene_id') %>%
    inner_join(e2g, by = 'gene_id') %>% 
    rename(EnsemblID = gene_id, GeneName = gene_name) %>%
    arrange(p.value) %>%
    select(EnsemblID, GeneName, AvgExpr, logFC, p.value, q.value)
  return(res)
}

# Agonist vs. Control
res_h1 <- output('Agonist')
fwrite(res_h1, './Results/Agonist_vs_Control.csv')
plot_md(res_h1, probes = 'EnsemblID',
        title = 'Mean-Difference Plot:\nAgonist vs. Control')
plot_volcano(res_h1, probes = 'EnsemblID',
             title = 'Volcano Plot:\nAgonist vs. Control')
top_genes <- res_h1$EnsemblID[seq_len(1000)]
tmp <- clin[Condition != 'Antagonist']
mat <- rld[rownames(rld) %in% top_genes, tmp$SampleID]
plot_heatmap(mat, group = list(Condition = tmp$Condition),
             title = 'Top Genes:\nAgonist vs. Control')

# Antagonist vs. Control
res_h2 <- output('Antagonist')
fwrite(res_h2, './Results/Antagonist_vs_Control.csv')
plot_md(res_h2, probes = 'EnsemblID',
        title = 'Mean-Difference Plot:\nAntagonist vs. Control')
plot_volcano(res_h2, probes = 'EnsemblID',
             title = 'Volcano Plot:\nAntagonist vs. Control')
top_genes <- res_h2$EnsemblID[seq_len(1000)]
tmp <- clin[Condition != 'Agonist']
mat <- rld[rownames(rld) %in% top_genes, tmp$SampleID]
plot_heatmap(mat, group = list(Condition = tmp$Condition),
             title = 'Top Genes:\nAntagonist vs. Control')










