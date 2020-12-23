
# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read in count table
read.counts <- read.table(file = "Filtered_sfnorm_counts.txt", header = TRUE)
names(read.counts)

# setwd("~/Documents/Laboratory/RNA_seq/ARSER-master/")

# get counts for WT2U
wt2u.counts <- cbind(row.names(read.counts), read.counts[ , c(37:48)])
head(wt2u.counts)
names(wt2u.counts) <- c("FLYBASE", seq(0, 44, by = 4))
head(wt2u.counts)
write.table(wt2u.counts, file = "ARSER_WT2U.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# get counts for WT2U_old
wt2uold.counts <- cbind(row.names(read.counts), read.counts[ , c(25:36)])
head(wt2uold.counts)
names(wt2uold.counts) <- c("FLYBASE", seq(0, 44, by = 4))
head(wt2uold.counts)
write.table(wt2uold.counts, file = "ARSER_WT2U_old.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# get counts for INC1
inc1.counts <- cbind(row.names(read.counts), read.counts[ , c(1:12)])
head(inc1.counts)
names(inc1.counts) <- c("FLYBASE", seq(0, 44, by = 4))
head(inc1.counts)
write.table(inc1.counts, file = "ARSER_INC1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# get counts for WAKED2
waked2.counts <- cbind(row.names(read.counts), read.counts[ , c(13:24)])
head(waked2.counts)
names(waked2.counts) <- c("FLYBASE", seq(0, 44, by = 4))
head(waked2.counts)
write.table(waked2.counts, file = "ARSER_WAKED2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## Hughes data-----------------------------------------------------------
# get counts for WT2U
# wt2u.counts <- cbind(row.names(read.counts), read.counts[ , c(1:8)])
# head(wt2u.counts)
# names(wt2u.counts) <- c("FLYBASE", seq(0, 42, by = 6))
# head(wt2u.counts)
# write.table(wt2u.counts, file = "ARSER_Hughes_CS.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# per.counts <- cbind(row.names(read.counts), read.counts[ , c(9:16)])
# head(per.counts)
# names(per.counts) <- c("FLYBASE", seq(0, 42, by = 6))
# head(per.counts)
# write.table(per.counts, file = "ARSER_Hughes_PER.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 

