library(scRNAseq)
sce <- BachMammaryData(location=FALSE)

# Saving in version 3 of the 10x format.
library(DropletUtils)
write10xCounts("mtx", assay(sce), paste0(sce$Sample, ".", sce$Barcode), gene.symbol=rowData(sce)$Symbol, version="3", overwrite=TRUE)

# Saving in version 3 of the 10x format.
library(HDF5Array)
name <- "tenx.h5"
unlink(name)
writeTENxMatrix(assay(sce), filepath=name, group="matrix")
rhdf5::h5createGroup(name, "matrix/features")
rhdf5::h5write(rownames(sce), name, "matrix/features/id")
rhdf5::h5write(rowData(sce)$Symbol, name, "matrix/features/name")

# Saving a H5AD file.
library(zellkonverter)
assay(sce) <- as(assay(sce), "dgCMatrix")
name <- "csc.h5ad"
unlink(name)
writeH5AD(sce, file=name)
