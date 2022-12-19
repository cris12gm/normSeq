args = commandArgs(trailingOnly=TRUE)
library(sva)
library(stringr)
library(RColorBrewer)
library(ggplot2)

matrix <- read.table(args[1],header=TRUE,sep="\t")

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

batchSheet <- read.table(args[2],header=T,sep="\t")
batchSheet$Sample <- make.names(batchSheet$sample)

row.names(batchSheet) <- batchSheet$sample
batchSheet$sample <- NULL

# output <- args[3]



# png(str_replace(output, ".png", "_single_row.png"),res=300,width=1800,height=1200)

# pheatmap(df, scale = "row", show_rownames=F, symm = T,
#         clustering_method="single",show_colnames = T,
#         cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(output,res=300,width=1800,height=1200)
# row_single <- pheatmap(df, scale = "row", show_rownames=F, symm = T,
#                         clustering_method="single",show_colnames = T,
#                         cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_complete_row.png"),res=300,width=1800,height=1200)
# row_complete <- pheatmap(df, scale = "row", show_rownames=F, symm = T,
#                        clustering_method="complete",show_colnames = T,
#                        cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_average_row.png"),res=300,width=1800,height=1200)
# row_average <- pheatmap(df, scale = "row", show_rownames=F, symm = T,
#                          clustering_method="average",show_colnames = T,
#                          cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_median_row.png"),res=300,width=1800,height=1200)
# row_median <- pheatmap(df, scale = "row", show_rownames=F, symm = T,
#                         clustering_method="median",show_colnames = T,
#                         cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_single_col.png"),res=300,width=1800,height=1200)
# col_single <- pheatmap(df, scale = "column", show_rownames=F, symm = T,
#                        clustering_method="single",show_colnames = T,
#                        cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_complete_col.png"),res=300,width=1800,height=1200)
# col_complete <- pheatmap(df, scale = "column", show_rownames=F, symm = T,
#                          clustering_method="complete",show_colnames = T,
#                          cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_average_col.png"),res=300,width=1800,height=1200)
# col_average <- pheatmap(df, scale = "column", show_rownames=F, symm = T,
#                         clustering_method="average",show_colnames = T,
#                         cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_median_col.png"),res=300,width=1800,height=1200)
# col_median <- pheatmap(df, scale = "column", show_rownames=F, symm = T,
#                        clustering_method="median",show_colnames = T,
#                        cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_single_none.png"),res=300,width=1800,height=1200)
# none_single <- pheatmap(df, scale = "none", show_rownames=F, symm = T,
#                        clustering_method="single",show_colnames = T,
#                        cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_complete_none.png"),res=300,width=1800,height=1200)
# none_complete <- pheatmap(df, scale = "none", show_rownames=F, symm = T,
#                           clustering_method="complete",show_colnames = T,
#                           cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_average_none.png"),res=300,width=1800,height=1200)
# none_average <- pheatmap(df, scale = "none", show_rownames=F, symm = T,
#                          clustering_method="average",show_colnames = T,
#                          cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()

# png(str_replace(output, ".png", "_median_none.png"),res=300,width=1800,height=1200)                        
# none_median <- pheatmap(df, scale = "none", show_rownames=F, symm = T,
#                         clustering_method="median",show_colnames = T,
#                         cluster_rows = F, cluster_cols = T,annotation_col=samplesheet)
# dev.off()