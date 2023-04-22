args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(RUVSeq))
suppressMessages(library(EDASeq))

matrix <- read.table(args[1],header=TRUE,sep="\t")

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

annotationSheet <- read.table(args[2],header=T,sep="\t",check.names=FALSE)
row.names(annotationSheet) <- annotationSheet$sample
annotationSheet$sample <- NULL


x<-as.factor(annotationSheet$replicate)
set<-newSeqExpressionSet(as.matrix(df), phenoData=data.frame(x,row.names=colnames(df)))

genes<-rownames(df)
differences<-makeGroups(x)

set2<-RUVs(set,genes,k=1,differences)

normalized <- as.data.frame(normCounts(set2))
outfile <- args[3]

normalized <- tibble::rownames_to_column(normalized, "name") # Apply rownames_to_column

write.table(normalized,outfile,sep="\t",row.names=FALSE, quote=F,col.names=TRUE)