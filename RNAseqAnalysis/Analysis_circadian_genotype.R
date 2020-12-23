
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

wt2u.union <- read.table(file = "circadian_union_JTK_WT2U_result.txt", header = TRUE)
inc1.union <- read.table(file = "circadian_union_JTK_INC1_result.txt", header = TRUE)
waked2.union <- read.table(file = "circadian_union_JTK_WAKED2_result.txt", header = TRUE)
wt2uold.union <- read.table(file = "circadian_union_JTK_WT2U_old_result.txt", header = TRUE)

head(wt2u.union)

wt2u.common <- intersect(row.names(wt2u.union), row.names(wt2uold.union))
wt2u.total <- union(row.names(wt2u.union), row.names(wt2uold.union))
wt2u.young <- setdiff(row.names(wt2u.union), row.names(wt2uold.union))
wt2u.old <- setdiff(row.names(wt2uold.union), row.names(wt2u.union))

inc1.common <- intersect(row.names(inc1.union), wt2u.common)
inc1.young <- intersect(row.names(inc1.union), wt2u.young)
inc1.old <- intersect(row.names(inc1.union), wt2u.old)
inc1.only <- setdiff(row.names(inc1.union), wt2u.total)

waked2.common <- intersect(row.names(waked2.union), wt2u.common)
waked2.young <- intersect(row.names(waked2.union), wt2u.young)
waked2.old <- intersect(row.names(waked2.union), wt2u.old)
waked2.only <- setdiff(row.names(waked2.union), wt2u.total)

mutant.common <- intersect(row.names(inc1.union), row.names(waked2.union))

rpkmPlotbyID("FBgn0036534")


intersect.all <- intersect(intersect(intersect(row.names(wt2u.union), row.names(wt2uold.union)), row.names(inc1.union)), row.names(waked2.union))

anno.DGE <- select(org.Dm.eg.db, keys = intersect.all, keytype = "FLYBASE", columns=c("FLYBASECG", "SYMBOL","GENENAME"))

union.all <- union(union(union(row.names(wt2u.union), row.names(wt2uold.union)), row.names(inc1.union)), row.names(waked2.union))

intersect.all <- read.counts[intersect.all,]
union.all <- read.counts[union.all,]

jtk.out <- read.table(file = "JTK_CYCLE_WT2U_out.txt", header = TRUE)
row.names(jtk.out) <- jtk.out$x


# extract subset of information from jtk.out
phase <- jtk.out[row.names(union.all), ]

# order phase and data according to phase
phase <- phase[order(phase$LAG),]
circadian.data <- union.all[row.names(phase),]

# draw heatmap for each genotype
for (i in c(1,13,25,37)){
  plot.data <- circadian.data[, i:(i+11)]
  aheatmap(plot.data, Rowv = NA, Colv = NA, scale = "row")
  # aheatmap(plot.data, Rowv = NA, Colv = NA, scale = "row", filename = paste("heatmap_all_circadian_", i, ".pdf", sep=""))
}

# save wt2u plot
wt2u.data <- circadian.data[, 37:48]
aheatmap(wt2u.data, Rowv = NA, Colv = NA, scale = "row", filename = "heatmap_all_circadian.pdf")

# save ordered circadian data with JTK result
out.df <- cbind(phase[ , 2:6], circadian.data)
#write.table(out.df, file = "circadian_union_JTK_WT2U.txt", sep = "\t", quote = FALSE, row.names = TRUE)
