args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))

matrix <- read.table(args[1],header=TRUE,sep="\t")

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

# dds <- DESeqDataSetFromMatrix(countData = df)

# print(dds)

# dds <- DESeq(dds)
# resultsNames(dds) # lists the coefficients
# res <- results(dds, name="condition_trt_vs_untrt")
# # or to shrink log fold changes association with condition:
# res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")



# #create DGEList object
# d <- DGEList(counts=df)
# d <- estimateCommonDisp(d)

# method <- args[2]

# if (method == "UQ") {
#     uq <- calcNormFactors(d, method="upperquartile")
#     normalized <- as.data.frame(uq$counts)
# }

# if (method == "TMM") {
#     tmm <- calcNormFactors(d, method="TMM")
#     normalized <- as.data.frame(tmm$counts)
# }
# if (method == "RLE") {
#     rle <- calcNormFactors(d, method="RLE")
#     normalized <- as.data.frame(rle$counts)
# }

# outfile <- args[3]
#                                          # Duplicate example data
# normalized <- tibble::rownames_to_column(normalized, "name") # Apply rownames_to_column

# write.table(normalized,outfile,sep="\t",row.names=FALSE, quote=F,col.names=TRUE)