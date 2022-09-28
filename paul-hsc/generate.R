library(scRNAseq)
sce <- PaulHSCData(location=FALSE)

# Saving in version 3 of the 10x format.
library(DropletUtils)
write10xCounts("mtx", assay(sce), version="3", overwrite=TRUE)

# Saving in version 3 of the 10x format.
library(HDF5Array)
source("../_scripts/compact10x.R")
name <- "tenx.h5"
unlink(name)
compact10x(name, assay(sce), rownames(sce), rownames(sce))

# Saving RDS files. We save this one as a base SE because
# there's not much going on anyway.
se <- SummarizedExperiment(list(counts=assay(sce)), colData=colData(sce), rowData=rowData(sce))
saveRDS(se, "se.rds")
