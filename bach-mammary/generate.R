library(scRNAseq)
sce <- BachMammaryData(location=FALSE)

# Saving in version 3 of the 10x format.
library(DropletUtils)
write10xCounts("mtx", assay(sce), paste0(sce$Sample, ".", sce$Barcode), gene.symbol=rowData(sce)$Symbol, version="3", overwrite=TRUE)

# Saving in version 3 of the 10x format.
source("../_scripts/compact10x.R")
name <- "tenx.h5"
unlink(name)
compact10x(name, assay(sce), id=rownames(sce), symbol=rowData(sce)$Symbol)
