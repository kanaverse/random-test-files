# This creates a nice small subset for tests to play with.
# It assumes that 'download.sh' has already been run.

library(DropletUtils)
se <- read10xCounts("immune_3.0.0-tenx.h5")

source("../_scripts/compact10x.R")

sub <- se[,1:3000]
path <- "immune_3.0.0_sub-tenx.h5"
unlink(path)
compact10x(path, assay(sub), rownames(sub), rowData(sub)$Symbol) 
h5write(rowData(sub)$Type, path, "matrix/features/feature_type")

path <- "immune_3.0.0_sub"
unlink(path, recursive=TRUE)
DropletUtils::write10xCounts(path, as(assay(sub), "dgCMatrix"), barcodes=colData(sub)$Barcode, gene.symbol=rowData(sub)$Symbol, gene.type=rowData(sub)$Type, version="3")

file.rename(file.path(path, "matrix.mtx.gz"), paste0(path, "-matrix.mtx.gz"))
file.rename(file.path(path, "barcodes.tsv.gz"), paste0(path, "-barcodes.tsv.gz"))
file.rename(file.path(path, "features.tsv.gz"), paste0(path, "-features.tsv.gz"))
