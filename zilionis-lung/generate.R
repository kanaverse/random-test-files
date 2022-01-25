library(scRNAseq)
sce <- ZilionisLungData()

# Saving in version 3 of the 10x format.
library(DropletUtils)
write10xCounts("mtx", assay(sce), version="3", overwrite=TRUE)

# Saving in version 3 of the 10x format.
source("../_scripts/compact10x.R")
name <- "tenx.h5"
unlink(name)
compact10x(name, assay(sce), rownames(sce), rownames(sce))
