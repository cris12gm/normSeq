args = commandArgs(trailingOnly=TRUE)

library(edgeR)
library(dplyr)

matrix <- read.table(args[1],header=TRUE,sep="\t",check.names=FALSE)

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
df <- matrix[keep, ]

method <- args[2]
annotation <- args[3]

pvalue <- args[4]
group1 <- args[5]
group2 <- args[6]

sampleSheet <- read.table(annotation,header=T,sep="\t",check.names=FALSE)
groups <- sampleSheet$group
sampletypevalues <- groups[!duplicated(groups)]  # Getting the group levels

# #create DGEList object
d <- DGEList(counts=df,group=factor(groups))
d <- estimateCommonDisp(d)

if (method == "UQ") {
    edgeR_table <- calcNormFactors(d, method="upperquartile")
}

if (method == "TMM") {
    edgeR_table <- calcNormFactors(d, method="TMM")
}
if (method == "RLE") {
    edgeR_table <- calcNormFactors(d, method="RLE")
}

edgeR_table <- estimateCommonDisp(edgeR_table)
  
edgeR_table <- estimateTagwiseDisp(edgeR_table)
selected <- NULL
dgeTest <- exactTest(edgeR_table, pair=c(group1,group2))

tt <- topTags(dgeTest, n=Inf)
selected_samples <- (which(groups==group1 | groups==group2))

group1Element <- (which(groups==group1))
group2Element <- (which(groups==group2))

# Obtaining the pseudocounts for each sample
z <- as.data.frame(edgeR_table$pseudo.counts)[ ,selected_samples]
  # Obtaining the statistical results from the edge test
t <- as.data.frame(tt)
  # Merging the above info into one table
data <- merge(z, t, by="row.names")
row.names(data) <- data$Row.names  # row names manipulation
data$Row.names <- NULL  # row names manipulation

group1_Mean = rowMeans(edgeR_table$pseudo.counts[,group1Element])
group2_Mean = rowMeans(edgeR_table$pseudo.counts[,group2Element])

data$group1Mean <- group1_Mean
data$group2Mean <- group2_Mean

names(data)[names(data) == 'group1Mean'] <- gsub(" ","",paste(group1,"_mean"))
names(data)[names(data) == 'group2Mean'] <- gsub(" ","",paste(group2,"_mean"))

# Selecting only features with FDR lower than the input p-value
selected <- which(data$FDR<=pvalue)
  # Obtaining the final matrix of selected features
result <- data[selected, ]

output <- args[7]
result <- tibble::rownames_to_column(as.data.frame(result), "name")

write.table(result,output,sep="\t",row.names=FALSE, quote=FALSE,col.names=TRUE)