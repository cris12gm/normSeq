args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(EDASeq))
suppressMessages(library(RColorBrewer))


matrix <- read.table(args[1],header=TRUE,sep="\t",check.names=FALSE,row.names=1)

annotationSheet <- read.table(args[2],header=T,sep="\t",check.names=FALSE)
row.names(annotationSheet) <- annotationSheet$sample
annotationSheet$sample <- NULL

set<-newSeqExpressionSet(as.matrix(matrix), phenoData=data.frame(as.factor(annotationSheet$group),row.names=colnames(matrix)))

#To get colors
num_uniq <- unique(annotationSheet$group)
maxColor <- length(num_uniq)

if (maxColor<=9){
  cols<-brewer.pal(maxColor,"Set1")
}else{
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  cols <- sample(col_vector, maxColor)
}


groups <- as.data.frame(cbind(annotationSheet$group))
groups$V1 <- as.factor(groups$V1)

colourGroup <- cols[groups$V1]
ann <- cbind(colourGroup)

output <- args[3]
png(output,res=300,width=1000,height=900)
par("mar" = c(5, 3, 1, 1),"mgp" = c(1.5, 1, 0))
plotRLE(set,outline=FALSE,col=ann,las=2,cex.axis=0.5,font="Palatino")+title(ylab="RLE",cex.lab=0.8)
dev.off()