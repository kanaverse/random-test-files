library(scRNAseq)
library(alabaster.sce)
alabaster.base::saveDataFrameFormat("hdf5")
sce <- ZeiselBrainData()

dump <- function(sce, path, ...) {
    unlink(path, recursive=TRUE)
    dir.create(path)
    info <- stageObject(sce, path, ".", ...)
    .writeMetadata(info, dir=path)

    tfile <- paste0(path, ".tar.gz")
    unlink(tfile)
    tar(tfile, files=list.files(path, full=TRUE, recursive=TRUE), compression="gzip")
}

local({
    assay(sce) <- as(assay(sce), "dgCMatrix")
    dump(sce, "zeisel-brain-sparse")

    mainExpName(sce) <- NULL 
    altExps(sce) <- list()
    dump(sce, "zeisel-brain-stripped")
})

local({
    assay(sce) <- as.matrix(assay(sce))
    dump(sce, "zeisel-brain-dense")
})

local({
    library(scrapper)
    res <- analyze(as(counts(sce), "dgCMatrix"), rna.subsets=list(Mito=grepl("^mt-", rownames(sce))))
    sce2 <- convertAnalyzeResults(res)
    dump(sce2, "zeisel-brain-sparse-results")

    res <- analyze(as.matrix(counts(sce)), adt.x=as.matrix(counts(altExp(sce))), rna.subsets=list(Mito=grepl("^mt-", rownames(sce))))
    sce2 <- convertAnalyzeResults(res)
    dump(sce2, "zeisel-brain-dense-multimodal-results")
})
