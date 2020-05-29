
# setup
# setwd("~/Documents/Laboratory/RNA_seq")
# library(org.Dm.eg.db)

analysis <- "DESeq2_TC_"

# read in all three files
inc1.list <- read.table(file = paste(analysis,"INC1",".txt",sep=""), header = TRUE)
waked2.list <-  read.table(file = paste(analysis,"WAKED2",".txt",sep=""), header = TRUE)
wt2uold.list <-  read.table(file = paste(analysis,"WT2U_old",".txt",sep=""), header = TRUE)

# get row names of each file
inc1.genename <- row.names(inc1.list)
waked2.genename <- row.names(waked2.list)
wt2uold.genename <- row.names(wt2uold.list)

intersect.genename <- intersect(intersect(inc1.genename, waked2.genename), wt2uold.genename)

# extract only the intersected data, modify column names
inc1.list <- inc1.list[intersect.genename,]
colnames(inc1.list) <- c(paste(colnames(inc1.list),"INC1",sep="_"))
waked2.list <- waked2.list[intersect.genename,]
colnames(waked2.list) <- c(paste(colnames(waked2.list),"WAKED2",sep="_"))
wt2uold.list <- wt2uold.list[intersect.genename,]
colnames(wt2uold.list) <- c(paste(colnames(wt2uold.list),"WT2U_old",sep="_"))

# write intersected data to file
candidate.list <- cbind(inc1.list, waked2.list, wt2uold.list)
write.table(candidate.list, file = paste(analysis,"intersect",".txt",sep=""), sep = "\t", quote = FALSE, row.names = TRUE)

# write intersected annotation to file
# candidates.annotation <- select(org.Dm.eg.db, keys = intersect.genename, keytype = "FLYBASE", columns=c("FLYBASECG", "GENENAME"))
# write.table(cbind(candidates.annotation,candidate.list), file = paste(analysis,"intersect_anno",".txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)


# plot all the genes when dseq.ds is already available
# for (genename in intersect.genename){
#   plotCounts(dseq.ds, gene = genename, normalized = TRUE)
# }



