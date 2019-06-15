# Set working directory
setwd('./Documents/QMUL/Crohns')

# Load libraries
library(AnnotationHub)
library(ensembldb)

# Load genome build 
ah <- AnnotationHub()
hs91 <- query(ah, c('ensembl', 'gtf', '91', 'Homo sapiens'))
gr <- hs91[[3]]

# Create ensembledb
db <- ensDbFromGRanges(gr, organism = 'Homo_sapiens',
                       genomeVersion = 'GRCh38', version = 91)
edb <- EnsDb(db)

# Create, export t2g df
t2g <- transcripts(edb, columns = c('tx_id', 'gene_id', 'gene_name'),
                   return.type = 'data.frame')
saveRDS(t2g, './Data/Hs91.t2g.rds')

# Now clean up the kallisto files
library(data.table)
clin <- fread('./Data/clinical.csv')
files <- file.path('./Data/Counts', clin$SampleID, 'abundance.tsv')
for (file in files) {
  tmp <- fread(file)
  tmp[, target_id := gsub('\\..*', '', target_id)]
  fwrite(tmp, file)
}
