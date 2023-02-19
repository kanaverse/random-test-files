library(scRNAseq)
sce <- ZeiselBrainData()
sub <- sce[1:1000,]
library(zellkonverter)
library(HDF5Array)
source("../_scripts/compact10x.R")

for (s in c(FALSE, TRUE)) {
    if (s) {
        x <- sub
        suffix <- "sub."
    } else {
        x <- sce
        suffix <- ""
    }

    # Saving in Matrix Market.
    library(DropletUtils)
    write10xCounts(paste0(suffix, "mtx"), as(assay(x),"dgCMatrix"), version="3", overwrite=TRUE)

    # Saving in dense format.
    name <- sprintf("dense.%sh5", suffix)
    unlink(name)
    writeHDF5Array(assay(x), filepath=name, name="matrix", chunkdim=c(1000, 500))

    # Saving in version 3 of the 10x format.
    name <- sprintf("tenx.%sh5", suffix)
    unlink(name)
    compact10x(name, assay(x), rownames(x), rownames(x))

    # Saving in H5AD format.
    assay(x) <- as(assay(x), "dgCMatrix")
    name <- sprintf("csc.%sh5ad", suffix)
    unlink(name)
    writeH5AD(x, file=name)

    # Saving in RDS format.
    if (!s) {
        saveRDS(x, sprintf("csc.%srds", suffix))
    }
}

# Saving as an RDS with GRanges.
sce2 <- ZeiselBrainData(ensembl=TRUE, location=TRUE)
saveRDS(sce2, sprintf("dense_plus_positions.rds"))

# Saving it after a little bit of an analysis.
library(scran.chan)
res <- quickBasicAnalysis(counts(sce), qc.subsets=list(Mito=grepl("^mt-", rownames(sce))))
sce3 <- marshalToSCE(sce, res)

counts(sce3) <- as(counts(sce3), "dgCMatrix") # to save some space
logcounts(sce3) <- as(logcounts(sce3), "dgCMatrix") # for simplicity for now.
saveRDS(sce3, "with_results.rds")

zellkonverter::writeH5AD(sce3, file="with_results.h5ad")
