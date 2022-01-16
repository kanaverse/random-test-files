library(scRNAseq)
sce <- ZeiselBrainData()

library(HDF5Array)
writeHDF5Array(assay(sce), filepath="zeisel.dense.h5", name="matrix", chunkdim=c(1000, 500))

# Saving in version 3 of the 10x format.
writeTENxMatrix(assay(sce), filepath="zeisel.tenx.h5", group="matrix")
rhdf5::h5createGroup("zeisel.tenx.h5", "matrix/features")
rhdf5::h5write(rownames(sce), "zeisel.tenx.h5", "matrix/features/id")
rhdf5::h5write(rownames(sce), "zeisel.tenx.h5", "matrix/features/name")

library(zellkonverter)
assay(sce) <- as(assay(sce), "dgCMatrix")
writeH5AD(sce, file="zeisel.csc.h5ad")


