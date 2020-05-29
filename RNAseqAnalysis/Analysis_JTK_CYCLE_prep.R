
# set working directory
setwd("~/Documents/Laboratory/RNA_seq")

# read in count table
read.counts <- read.table(file = "Filtered_sfnorm_counts.txt", header = TRUE)
names(read.counts)

# get the first column
gene.list <- row.names(read.counts)
head(gene.list)
# note this command puts an "x" on top of the file
write.table(gene.list, file = "JTK_CYCLE_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# get counts for WT2U
wt2u.counts <- read.counts[ , c(37:48)]
head(wt2u.counts)
write.table(wt2u.counts, file = "JTK_CYCLE_WT2U.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# get counts for WT2U_old
wt2uold.counts <- read.counts[ , c(25:36)]
head(wt2uold.counts)
write.table(wt2uold.counts, file = "JTK_CYCLE_WT2U_old.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# get counts for INC1
inc1.counts <- read.counts[ , c(1:12)]
head(inc1.counts)
write.table(inc1.counts, file = "JTK_CYCLE_INC1.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# get counts for WAKED2
waked2.counts <- read.counts[ , c(13:24)]
head(waked2.counts)
write.table(waked2.counts, file = "JTK_CYCLE_WAKED2.txt", sep = "\t", quote = FALSE, row.names = TRUE)


