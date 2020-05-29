
# set up
setwd("~/Documents/Laboratory/RNA_seq")

# read and combind feature counts tables
pool1.counts <- read.table(file = "counts_pool1.matrix", header = TRUE)
pool2.counts <- read.table(file = "counts_pool2.matrix", header = TRUE)
pool3.counts <- read.table(file = "counts_pool3.matrix", header = TRUE)

# set Geneid as row names
row.names(pool1.counts) <- pool1.counts$Geneid
row.names(pool2.counts) <- pool2.counts$Geneid
row.names(pool3.counts) <- pool3.counts$Geneid

# get genome information for each gene
counts.info <- pool1.counts[ , c(2:6)]

# save only count numbers
pool1.counts <- pool1.counts[ , c(7:22)]
pool2.counts <- pool2.counts[ , c(7:22)]
pool3.counts <- pool3.counts[ , c(7:22)]

# put data from three files into one
read.counts <- cbind(pool1.counts, pool2.counts, pool3.counts)
read.counts <- read.counts[order(names(read.counts))]
names(read.counts)

# set the appropriate column names
names(read.counts) <- c(c(paste("INC1",c("00", "12", "16", "20", "24", "28", "32", "36", "04", "40", "44", "08"),sep="_ZT")), c(paste("WAKED2",c("00", "12", "16", "20", "24", "28", "32", "36", "04", "40", "44", "08"),sep="_ZT")), c(paste("WT2U_old",c("00", "12", "16", "20", "24", "28", "32", "36", "04", "40", "44", "08"),sep="_ZT")), c(paste("WT2U",c("00", "12", "16", "20", "24", "28", "32", "36", "04", "40", "44", "08"),sep="_ZT")))
read.counts <- read.counts[order(names(read.counts))]
names(read.counts)

# write the counts and genome information to file
write.table(read.counts, file = "Raw_data_combined.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(counts.info, file = "Raw_data_combined_info.txt", sep = "\t", quote = FALSE, row.names = TRUE)
