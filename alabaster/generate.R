library(scRNAseq)
library(alabaster.sce)
sce <- ZeiselBrainData()

dump <- function(sce, path, ...) {
    unlink(path, recursive=TRUE)
    saveObject(sce, path, ...)
    tfile <- paste0(path, ".zip")
    unlink(tfile)
    zip(tfile, files=list.files(path, full=TRUE, recursive=TRUE))
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
    dump(sce2, "zeisel-brain-sparse-results-delayed", DelayedArray.preserve.ops=TRUE)
    dump(sce2, "zeisel-brain-sparse-results-delayed-external", DelayedArray.preserve.ops=TRUE, DelayedArray.force.external=TRUE)

    res <- analyze(as.matrix(counts(sce)), adt.x=as.matrix(counts(altExp(sce))), rna.subsets=list(Mito=grepl("^mt-", rownames(sce))))
    sce2 <- convertAnalyzeResults(res)
    dump(sce2, "zeisel-brain-dense-multimodal-results")
    dump(sce2, "zeisel-brain-dense-multimodal-results-delayed", DelayedArray.preserve.ops=TRUE)
    dump(sce2, "zeisel-brain-dense-multimodal-results-delayed-external", DelayedArray.preserve.ops=TRUE, DelayedArray.force.external=TRUE)
})
