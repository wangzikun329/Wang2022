
circadian.data <- read.table(file = "circadian_union_JTK_WT2U_result.txt", header = TRUE)
tc.data <- read.table(file = "DGE_TC_FINAL.txt", header = TRUE)
tc.direction <- read.csv(file = "DGE_TC_FINAL_same_direction.csv", header = TRUE)
tc.fc <- read.csv(file = "DGE_TC_FINAL_same_direction_min_fc.csv", header = TRUE)

row.names(tc.direction) <- tc.direction$FLYBASE
row.names(tc.fc) <- tc.fc$FLYBASE

circadian.tc <- intersect(row.names(circadian.data), row.names(tc.data))
circadian.dir <- intersect(row.names(circadian.data), row.names(tc.direction))
circadian.fc <- intersect(row.names(circadian.data), row.names(tc.fc))

out.df <- circadian.data[circadian.tc, ]

# write.table(out.df, file = "DGE_circadian_intersect.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# order phase and data according to phase
circadian.tc <- out.df[order(out.df$LAG),]

# subset of DGE genes with same change direction
circadian.dir <- circadian.data[circadian.dir, ]
dir.circadian <- tc.data[row.names(circadian.dir),]
dir.circadian.up <- dir.circadian[dir.circadian$log2FoldChange_INC1>0,]
dir.circadian.dn <- dir.circadian[dir.circadian$log2FoldChange_INC1<0,]
circadian.dir.up <- circadian.dir[row.names(dir.circadian.up),]
circadian.dir.dn <- circadian.dir[row.names(dir.circadian.dn),]
circadian.dir.up <- circadian.dir.up[order(circadian.dir.up$LAG),]
circadian.dir.dn <- circadian.dir.dn[order(circadian.dir.dn$LAG),]
# draw heatmap for each genotype
# genotype <- c(rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12), rep("WT2U", 12))
# for (i in c(1,13,25,37)){
#   plot.data <- circadian.data[, i:(i+11)]
#   aheatmap(plot.data, Rowv = NA, Colv = NA, scale = "row")
#   # aheatmap(plot.data, Rowv = NA, Colv = NA, scale = "row", filename = paste("heatmap_all_circadian_", genotype[i], ".pdf", sep=""))
# }

# save wt2u plot
wt2u.data <- circadian.tc[, 42:53]
aheatmap(wt2u.data, Rowv = NA, Colv = NA, scale = "row")
aheatmap(wt2u.data, Rowv = NA, Colv = NA, scale = "row", filename = "heatmap_DGE_circadian_intersect.pdf")

gene.data <- circadian.tc[, 6:53]
aheatmap(gene.data, Rowv = NA, Colv = NA, scale = "row")
aheatmap(gene.data, Rowv = NA, Colv = NA, scale = "row", filename = "heatmap_DGE_circadian_intersect_allgroup.pdf")

gene.data <- rbind(circadian.dir.dn, circadian.dir.up)[, c(c(42:53), c(6:41))]
row.anno <- c(rep(0, nrow(circadian.dir.dn)), rep(1, nrow(circadian.dir.up)))
aheatmap(gene.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 5, cellwidth = 10,
         annRow = row.anno, annColors = "-RdBu")
aheatmap(gene.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 10, cellwidth = 10,
         annRow = row.anno, annColors = "-RdBu", filename = "heatmap_DGE__circadian_intersect_direction.pdf")
