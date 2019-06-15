# Set working directory
setwd('~/Documents/QMUL/Crohns')

# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(tidyverse)
library(bioplotr)

# Import data
clin <- fread('./Data/clinical.csv')
t2g <- readRDS('./Data/Hs91.t2g.rds')
anno <- fread('./Data/Hs.anno.csv')
files <- file.path('./Data/Kallisto', clin$SampleID, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)

# Make dds object
dds <- DESeqDataSetFromTximport(txi, colData = clin, 
                                design = ~ Patient + Condition)
saveRDS(dds, './Data/dds_full.rds')
idx <- rownames(dds)

# Normalise, filter, transform counts
dds <- estimateSizeFactors(dds)
mat <- counts(dds, normalized = TRUE)
keep <- rowSums(mat > 1) >= 10
rld <- assay(rlog(dds[keep, ]))
colnames(rld) <- clin$SampleID
saveRDS(rld, './Data/rld_full.rds')

# Dispersion
plot_dispersion(dds[keep, ])

# Density
plot_density(rld, group = list(Condition = clin$Condition), 
             xlab = 'rlog Transformed Counts')

# Sample similarity
plot_similarity(rld, group = list(Condition = clin$Condition, Sex = clin$Sex), 
                covar = list(Age = clin$Age, RIN = clin$RIN),
                pal_covar = c('viridis', 'Spectral'))

# PCA 
plot_pca(rld, group = list(Patient = clin$Patient), label = TRUE)

# Rerun without outliers
outliers <- c('1911B-D', '2502A-F', '2502A-CH')
clin <- fread('./Data/clinical.csv')
clin <- clin[!SampleID %in% outliers]
files <- file.path('./Data/Kallisto', clin$SampleID, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)
dds <- DESeqDataSetFromTximport(txi, colData = clin, 
                                design = ~ Patient + Condition)
saveRDS(dds, './Data/dds.rds')
idx <- rownames(dds)
dds <- estimateSizeFactors(dds)
mat <- counts(dds, normalized = TRUE)
keep <- rowSums(mat > 1) >= 10
rld <- assay(rlog(dds[keep, ]))
colnames(rld) <- clin$SampleID
saveRDS(rld, './Data/rld.rds')

# Quick load
outliers <- c('1911B-D', '2502A-F', '2502A-CH')
clin <- fread('./Data/clinical.csv')
clin <- clin[!SampleID %in% outliers]
dds <- readRDS('./Data/dds.rds')
idx <- rownames(dds)
rld <- readRDS('./Data/rld.rds')
anno <- fread('./Data/Hs.anno.csv')

# PCA 
plot_pca(rld, group = list(Patient = clin$Patient), label = TRUE)

# KPCA
plot_kpca(rld, group = list(Patient = clin$Patient), label = TRUE)

# t-SNE
plot_tsne(rld, group = list(Patient = clin$Patient), label = TRUE)

# Drivers
tmp <- clin %>% select(SampleID, Patient, Condition, Sex, Age, RIN)
plot_drivers(rld, tmp, index = 'SampleID', block = 'Patient', unblock = 'Age', 
             alpha = 0.05, p.adj = 'fdr', n.pc = 10)




















