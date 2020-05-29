
# read in DESeq2 results and edgeR results
results.DESeq2 <- read.table(file = "DESeq2_TC_intersect.txt", header = TRUE)
results.edgeR.qCML <- read.table(file = "edgeR_TC_intersect.txt", header = TRUE)
# results.edgeR.glm <- read.table(file = "edgeRresults_glm_intersection.txt", header = TRUE)

# find intersection in two result tables
# intersect.genename <- intersect(intersect(row.names(results.DESeq2), row.names(results.edgeR.qCML)), row.names(results.edgeR.glm))
intersect.genename <- intersect(row.names(results.DESeq2), row.names(results.edgeR.qCML))

# subset the results table
results.DESeq2 <- results.DESeq2[intersect.genename, ]
results.edgeR.qCML <- results.edgeR.qCML[intersect.genename, ]
# results.edgeR.glm <- results.edgeR.glm[intersect.genename, ]

# annotate final result
library(org.Dm.eg.db)
# keytypes(org.Dm.eg.db)
# make a batch retrieval for all DE genes
intersect.genename.annotation <- select(org.Dm.eg.db, keys = intersect.genename, keytype = "FLYBASE", columns=c("FLYBASECG", "GENENAME"))
row.names(intersect.genename.annotation) <- intersect.genename.annotation$FLYBASE

# merge two tables into one
results.intersect <- cbind(results.DESeq2, results.edgeR.qCML)
results.anno <- cbind(intersect.genename.annotation, results.DESeq2, results.edgeR.qCML)

write.table(results.intersect, file = "DGE_TC_FINAL.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(results.anno, file = "DGE_TC_FINAL_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
