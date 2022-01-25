library(scRNAseq)
sce <- ZeiselBrainData()
sub <- sce[1:1000,]
library(zellkonverter)
library(HDF5Array)

for (s in c(FALSE, TRUE)) {
    if (s) {
        x <- sub
        suffix <- "sub."
    } else {
        x <- sce
        suffix <- ""
    }

    name <- sprintf("zeisel.dense.%sh5", suffix)
    unlink(name)
    writeHDF5Array(assay(x), filepath=name, name="matrix", chunkdim=c(1000, 500))

    # Saving in version 3 of the 10x format.
    name <- sprintf("zeisel.tenx.%sh5", suffix)
    unlink(name)
    writeTENxMatrix(assay(x), filepath=name, group="matrix")
    rhdf5::h5createGroup(name, "matrix/features")
    rhdf5::h5write(rownames(x), name, "matrix/features/id")
    rhdf5::h5write(rownames(x), name, "matrix/features/name")

    assay(x) <- as(assay(x), "dgCMatrix")
    name <- sprintf("zeisel.csc.%sh5ad", suffix)
    unlink(name)
    writeH5AD(x, file=name)
}
