args = commandArgs(trailingOnly=TRUE)

library(DESeq2)
library(dplyr)


matrix <- read.table(args[1],header=TRUE,sep="\t")

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
matfile <- matrix[keep, ]

annotation <- args[2]

pvalue <- args[3]
min_t <- args[4]
group1 <- args[5]
group2 <- args[6]

sampleSheet <- read.table(annotation,header=T,sep="\t")
groups <- sampleSheet$group
sampletypevalues <- groups[!duplicated(groups)]  # Getting the group levels

# Designing the data's factors which indicate the experimental group for each sample
samplefactors <- data.frame(row.names=colnames(matfile), condition = factor(groups, levels=sampletypevalues))

# Constructing the data set object that DESeq2 needs to perform the analysis
dds = DESeqDataSetFromMatrix(countData = round(matfile), colData = samplefactors, design = ~condition)


# Calling DESeq method
dds <- DESeq(dds)
# Obtaining the normalised matrix
norvalues <- counts(dds, normalized=TRUE)
# Mean expression of the normalised data
meanexpression <- rowMeans(norvalues)
# Discarding features with mean expression lower than
# the user-given threshold 'min_t'
rcounts <- matfile[which(meanexpression >= min_t), ]
# Reconstructing the data set object that DESeq2 needs to perform the analysis
dds_filtered = DESeqDataSetFromMatrix(countData = round(rcounts), colData = samplefactors, design = ~condition)
# Calling DESeq method
dds_filtered <- DESeq(dds_filtered)

# Using the estimates the size factors to the filtered data
sizeFactors(dds_filtered) <- sizeFactors(dds)
# Normalisation of the filtered data
ncounts <- counts(dds_filtered, normalized=TRUE)

# Extracting the table containing the test results for each feature in the original count table
curresult <- results(dds_filtered)
selected <- NULL


# Selecting the samples
selected_samples <- (which(groups==sampletypevalues[1] | groups==sampletypevalues[2]))

group1Element <- (which(groups==group1))
group2Element <- (which(groups==group2))

# Obtaining the list of features with absolute log2FoldChange greater than 1 and p adjusted value lower than the user-input value
selected <- row.names(curresult)[(abs(curresult$log2FoldChange)>=1) & (curresult$padj<=pvalue)]
# Removing rows with missing values on columns
selected <- na.omit(selected)

# Combing the normalised data along with statistical analysis results ("log2FoldChange", "pval", "padj")
ncounts_selected <- cbind(ncounts[ , selected_samples], curresult$log2FoldChange, curresult$pvalue, curresult$padj)

# Naming the new columns
colnames(ncounts_selected) <- c(head(colnames(ncounts_selected), n=-3), "logFC", "PValue", "FDR")
# Obtaining the final matrix of selected features
result <- as.data.frame(ncounts_selected[selected, ])

try(result$group1Mean <- rowMeans(subset(result, select = group1Element)))
try(result$group2Mean <- rowMeans(subset(result, select = group2Element)))

names(result)[names(result) == 'group1Mean'] <- gsub(" ","",paste(group1,"_mean"))
names(result)[names(result) == 'group2Mean'] <- gsub(" ","",paste(group2,"_mean"))


output <- args[7]
result <- tibble::rownames_to_column(as.data.frame(result), "name")
write.table(result,output,sep="\t",row.names=FALSE, quote=FALSE,col.names=TRUE)