args = commandArgs(trailingOnly=TRUE)

library(RUVSeq)

matrix <- read.table(args[1],header=TRUE,sep="\t")

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

design <- model.matrix(~x, data=df)
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")


normalized <- normalizeQuantiles(df, ties=TRUE)

outfile <- args[2]

normalized <- tibble::rownames_to_column(normalized, "name") # Apply rownames_to_column

write.table(normalized,outfile,sep="\t",row.names=FALSE, quote=F,col.names=TRUE)