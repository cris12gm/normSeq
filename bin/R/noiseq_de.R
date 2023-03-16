args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(NOISeq))
suppressMessages(library(dplyr))


matrix <- read.table(args[1],header=TRUE,sep="\t", check.names=FALSE)

#Remove empty and 0s

matrix <- na.omit(matrix)
row.names(matrix) <- matrix$name
matrix = subset(matrix, select = -c(name) )

keep <- rowSums(matrix>0) > 0
matfile <- matrix[keep, ]

method <- args[2]
annotation <- args[3]

pvalue <- args[4]
min_t <- args[5]
group1 <- args[6]
group2 <- args[7]

sampleSheet <- read.table(annotation,header=T,sep="\t")
sampleSheet <- sampleSheet[sampleSheet$group==group1 | sampleSheet$group==group2,]

groups <- sampleSheet$group
samplesAnalysis <- sampleSheet$sample
sampletypevalues <- groups[!duplicated(groups)]  # Getting the group levels

#Filter df
matfile <- t(t(matfile)[samplesAnalysis,])

# Designing the data's factors which indicate the experimental group for each sample
samplefactors <- data.frame(row.names=colnames(matfile), condition = factor(groups, levels=sampletypevalues))
# Importing all necessary information into a NOISeq object
countdata <- readData(data = matfile, factors = samplefactors)
# Trimmed Mean of M values (TMM) normalisation to correct the sequencing depth bias
# No length is taken into account
TMMvalues = tmm(assayData(countdata)$exprs, long = 1000, lc = 0)
# Mean expression of the normalised data
meanexpression <- rowMeans(TMMvalues)
# Discarding features with mean expression lower than
# the user-given threshold 'min_t'
TMM_filtered <- matfile[which(meanexpression>=min_t), ]
# Recollecting the filtered data for the NOISeq analysis
countdata_filtered = readData(data = TMM_filtered, factors = samplefactors)
# TMM normalisation of the filtered data. No length is taken into account
TMMvalues_filtered = tmm(assayData(countdata_filtered)$exprs, long = 1000, lc = 0)
# Computing the differential expression between experimental conditions from the filtered read count data
NOISeq = noiseq(countdata_filtered, k=0.5, lc=0, norm="tmm", factor="condition")
qvalue = 1-as.double(pvalue)
mynoiseq.deg = degenes(NOISeq, q = qvalue, M = NULL)

# Extract the table containing the test results for each feature of the original count table
curresult <- mynoiseq.deg

# Selecting the samples
selected_samples <- (which(groups==sampletypevalues[1] | groups==sampletypevalues[2]))
group1Element <- (which(groups==group1))
group2Element <- (which(groups==group2))

# Obtaining the list of features with absolute log2FoldChange greater than 1 and p value lower than the user-input value
selected <- curresult[(abs(curresult$M)>=1),]
result <- subset(selected, select = -c(D,ranking) )
result$pval <- 1 - result$prob

result <- result[(result$pval)<=pvalue,]

output <- args[8]
result <- tibble::rownames_to_column(as.data.frame(result), "name")
colnames(result) <- c(head(colnames(result), n=-3), "logFC", "prob", "PValue")


write.table(result,output,sep="\t",row.names=FALSE, quote=FALSE,col.names=TRUE)