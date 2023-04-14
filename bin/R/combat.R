args = commandArgs(trailingOnly=TRUE)
library(sva)
library(stringr)
library(RColorBrewer)
library(ggplot2)

matrix <- read.table(args[1],header=TRUE,sep="\t",check.names=FALSE)

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

batchSheet <- read.table(args[2],header=T,sep="\t",check.names=FALSE)
batchSheet$Sample <- make.names(batchSheet$sample)
row.names(batchSheet) <- batchSheet$sample
batchSheet$sample <- NULL

annotationSheet <- read.table(args[3],header=T,sep="\t",check.names=FALSE)
row.names(annotationSheet) <- annotationSheet$sample
annotationSheet$sample <- NULL

merged <- merge(batchSheet, annotationSheet, by.x = 0, by.y = 0)
row.names(merged) <- merged$Row.names
merged$Row.names <- NULL

batchSheet <- merged[order(merged$Sample),]
df <- t(t(df)[order(row.names(t(df))), ])

output <- args[4]

batch <- batchSheet$batchEffect
groups <- batchSheet$group

adjusted <- ComBat_seq(df, batch=batch, group=groups)
adjusted <- tibble::rownames_to_column(as.data.frame(adjusted), "name")

write.table(adjusted,output,sep="\t",row.names=FALSE, quote=FALSE,col.names=TRUE)