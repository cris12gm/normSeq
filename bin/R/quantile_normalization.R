args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(limma))

matrix <- read.table(args[1],header=TRUE,sep="\t",check.names=FALSE)

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

normalized <- normalizeQuantiles(df, ties=TRUE)

outfile <- args[2]

normalized <- tibble::rownames_to_column(normalized, "name") # Apply rownames_to_column

write.table(normalized,outfile,sep="\t",row.names=FALSE, quote=F,col.names=TRUE)