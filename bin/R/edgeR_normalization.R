args = commandArgs(trailingOnly=TRUE)

library(edgeR)
library(dplyr)

matrix <- read.table(args[1],header=TRUE,sep="\t")

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

#create DGEList object
d <- DGEList(counts=df)

method <- args[2]

if (method == "UQ") {
    d <- calcNormFactors(d, method="upperquartile")
    d <- estimateCommonDisp(d)
    uq <- cpm(d, normalized.lib.sizes=TRUE, log=FALSE)
    normalized <- as.data.frame(uq)
}

if (method == "TMM") {
    d <- calcNormFactors(d, method="TMM")
    d <- estimateCommonDisp(d)
    tmm <- cpm(d, normalized.lib.sizes=TRUE, log=FALSE)
    normalized <- as.data.frame(tmm)
}
if (method == "RLE") {
    d <- calcNormFactors(d, method="RLE")
    d <- estimateCommonDisp(d)
    rle <- cpm(d, normalized.lib.sizes=TRUE, log=FALSE)
    normalized <- as.data.frame(rle)
}

outfile <- args[3]
                                         # Duplicate example data
normalized <- tibble::rownames_to_column(normalized, "name") # Apply rownames_to_column

write.table(normalized,outfile,sep="\t",row.names=FALSE, quote=F,col.names=TRUE)