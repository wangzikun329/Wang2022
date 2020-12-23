
# setup
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggplot2)
library(DESeq2)
library(org.Dm.eg.db)
library(GO.db)

# install required packages for R 3.6
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("org.Dm.eg.db")
# BiocManager::install("vsn")
# BiocManager::install("NMF")
# BiocManager::install("GO.db")
# BiocManager::install("edgeR")
# BiocManager::install("data.table")


# read in raw data  
read.counts <- read.table(file = "Filtered_raw_counts.txt", header = TRUE)
read.rpkm <- read.table(file = "Raw_data_combined_rpkm.txt", header = TRUE)

# define sample metadata
sample.info <- data.frame(condition = c(rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12), rep("WT2U", 12)), row.names = names(read.counts) )
ZT <- rep(c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"), 4)
genotype = c(rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12), rep("WT2U", 12))

# generate DESeqDataSet
dseq.plot <- DESeqDataSetFromMatrix(count = read.counts, colData = sample.info, design = ~ condition)

# run DESeq
dseq.plot <- DESeq(dseq.plot)

# single gene plot with zt as x
retrieveName <- function(FBgnID){
  genename <- select(org.Dm.eg.db, keys = c(FBgnID, "FBgn0003996"), keytype = "FLYBASE", columns=c("GENENAME"))
  # genename <- strsplit(genename[,2]," ")[[1]][1]
  return(genename[1,2])
}

countPlotbyID <- function(FBgnID, returnData=FALSE, savePlottoPDF=FALSE){
  gene.data <- plotCounts(dseq.plot, gene = FBgnID, normalized = TRUE, returnData=TRUE)
  gene.plot <- ggplot(gene.data, aes(x=ZT, y=count, color=condition, group=condition)) + 
    geom_point() + geom_line() + theme_linedraw() + ylab("Normalized count") + theme(panel.grid = element_blank()) +
    ggtitle(paste(FBgnID,retrieveName(FBgnID),sep="\n")) + theme(plot.title = element_text(hjust = 0.5))
  if(savePlottoPDF == TRUE){
    pdf(paste(FBgnID, ".pdf", sep = ""))
    print(gene.plot)
    dev.off()
  }
  print(gene.plot)
  if(returnData == TRUE){
    colnames(gene.data)[1] <- paste0(FBgnID, "_count")
    return(gene.data[,1,drop=FALSE])
  }
}

countPlotbyName <- function(genename){
  FBgnID <- select(org.Dm.eg.db, keys = genename, keytype = "GENENAME", columns=c("FLYBASE"))
  countPlotbyID(FBgnID[,2])
}

rpkmPlotbyID <- function(FBgnID, returnData=FALSE, savePlottoPDF=FALSE){
  gene.data <- as.data.frame(t(read.rpkm[FBgnID, ]))
  gene.plot <- ggplot(gene.data, aes(x=ZT, y=gene.data[, 1], color=genotype, group=genotype)) +
    geom_point() + geom_line() + theme_linedraw() + ylab("RPKM") + theme(panel.grid = element_blank()) +
    ggtitle(paste(FBgnID,retrieveName(FBgnID),sep="\n")) + theme(plot.title = element_text(hjust = 0.5))
  if(savePlottoPDF == TRUE){
    pdf(paste(FBgnID, ".pdf", sep = ""))
    print(gene.plot)
    dev.off()
  }
  print(gene.plot)
  if(returnData == TRUE){
    colnames(gene.data) <- paste0(FBgnID, "_RPKM")
    return(gene.data)
  }
}

# rpkmPlotbyID <- function(FBgnID, returnData=FALSE, savePlottoPDF=FALSE){
#   gene.data <- as.data.frame(t(read.rpkm[FBgnID, ]))
#   gene.plot <- ggplot(gene.data, aes(x=ZT, y=gene.data[, 1], color=genotype, group=genotype)) + 
#     geom_point() + geom_line() + theme_linedraw() + ylab("RPKM") + 
#     ggtitle(paste(FBgnID,retrieveName(FBgnID),sep="\n")) +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       panel.grid = element_blank(),
#       panel.border = element_blank(),
#       axis.line = element_line(size = 0.5),
#       axis.ticks = element_line(size = 0.5)
#     )
#   if(savePlottoPDF == TRUE){
#     pdf(paste(FBgnID, ".pdf", sep = ""))
#     print(gene.plot)
#     dev.off()
#   }
#   print(gene.plot)
#   if(returnData == TRUE){
#     colnames(gene.data) <- paste0(FBgnID, "_RPKM")
#     return(gene.data)
#   }
# }

rpkmPlotbyName <- function(genename, returnData=FALSE, savePlottoPDF=FALSE){
  FBgnID <- select(org.Dm.eg.db, keys = genename, keytype = "GENENAME", columns=c("FLYBASE"))
  rpkmPlotbyID(FBgnID[,2], returnData, savePlottoPDF)
}

rpkmPlotbyCG <- function(geneCG, returnData=FALSE, savePlottoPDF=FALSE){
  FBgnID <- select(org.Dm.eg.db, keys = geneCG, keytype = "FLYBASECG", columns=c("FLYBASE"))
  rpkmPlotbyID(FBgnID[,2], returnData, savePlottoPDF)
}

countPlotbyID("FBgn0036585")
countPlotbyName("white")

rpkmPlotbyID("FBgn0036585")
# rpkmPlotbyID("FBgn0036585", savePlottoPDF = TRUE)
rpkmPlotbyName("period", returnData = TRUE)

# for (genename in genename.list){
#   gene.data <- plotCounts(dseq.ds, gene = genename, normalized = TRUE, returnData=TRUE)
#   print(ggplot(gene.data, aes(x=ZT, y=count, color=condition, group=condition)) + geom_point() + geom_line() + theme_linedraw())
#   # plotCounts(dseq.ds, gene = genename, normalized = TRUE)
# }
