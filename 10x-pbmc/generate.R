library(TENxPBMCData)
sce3 <- TENxPBMCData("pbmc3k")
colnames(sce3) <- sce3$Barcode
assay(sce3) <- as(assay(sce3), "dgCMatrix")

sce4 <- TENxPBMCData("pbmc4k")
colnames(sce4) <- sce4$Barcode
assay(sce4) <- as(assay(sce4), "dgCMatrix")

# Saving in version 3 of the 10x format.
library(DropletUtils)
write10xCounts("mtx-pbmc3k", assay(sce3), gene.symbol=rowData(sce3)$Symbol, version="3", overwrite=TRUE)
write10xCounts("mtx-pbmc4k", assay(sce4), gene.symbol=rowData(sce4)$Symbol, version="3", overwrite=TRUE)

# Saving in version 3 of the 10x format.
source("../_scripts/compact10x.R")

name <- "pbmc3k-tenx.h5"
unlink(name)
compact10x(name, assay(sce3), rownames(sce3), rowData(sce3)$Symbol)

name <- "pbmc4k-tenx.h5"
unlink(name)
compact10x(name, assay(sce4), rownames(sce4), rowData(sce4)$Symbol)

# Combining into a single object.
common <- intersect(rownames(sce3), rownames(sce4))
x <- cbind(
    assay(sce3)[common,,drop=FALSE],
    assay(sce4)[common,,drop=FALSE]
)
write10xCounts("combined", x, 
    gene.symbol=rowData(sce3)[match(common, rownames(sce3)),"Symbol"],
    version="3", overwrite=TRUE)

Y <- read.delim('combined/barcodes.tsv.gz', header=FALSE)
Y$V2 <- rep(c("3k", "4k"), c(ncol(sce3), ncol(sce4)))
gzhandle <- gzfile("combined/barcodes.tsv.gz", open="wb")
write.table(gzhandle, x=Y, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
close(gzhandle)
