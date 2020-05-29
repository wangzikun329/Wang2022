
genotype <- "TC_INC1"

# read in DESeq2 results and edgeR results
results.DESeq2 <- read.table(file = paste("DESeq2_",genotype,".txt",sep=""), header = TRUE)
results.edgeR <- read.table(file = paste("edgeR_",genotype,".txt",sep=""), header = TRUE)

# find intersection in two result tables
intersect.genename <- intersect(row.names(results.DESeq2), row.names(results.edgeR))

# subset the results table
results.DESeq2 <- results.DESeq2[intersect.genename, ]
results.edgeR <- results.edgeR[intersect.genename, ]

# annotate final result
# library(org.Dm.eg.db)
# keytypes(org.Dm.eg.db)
# make a batch retrieval for all DE genes
# intersect.genename.annotation <- select(org.Dm.eg.db, keys = intersect.genename, keytype = "FLYBASE", columns=c("FLYBASECG", "GENENAME"))
# row.names(intersect.genename.annotation) <- intersect.genename.annotation$FLYBASE

# merge two tables into one
results.intersect <- cbind(results.DESeq2, results.edgeR)
# results.anno <- cbind(intersect.genename.annotation, results.DESeq2, results.edgeR)

write.table(results.intersect, file = paste("DGE_",genotype,".txt",sep=""), sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(results.anno, file = paste("DGE_",genotype,"_anno.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
