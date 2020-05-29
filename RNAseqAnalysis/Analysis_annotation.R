
# the output file of this script could not be read by R
file.name <- "circadian_setdiff_RhythmicityLoss.txt"

# setup
setwd("~/Documents/Laboratory/RNA_seq")
library(org.Dm.eg.db)

# read in file with Flybase ID as row names
read.file <- read.table(file = file.name, header = TRUE)

# retrieve annotation from database
anno.DGE <- select(org.Dm.eg.db, keys = row.names(read.file), keytype = "FLYBASE", columns=c("FLYBASECG", "SYMBOL","GENENAME"))
out.df <- cbind(anno.DGE, read.file)

write.table(out.df, file = paste(file.name,"anno.txt", sep="."), sep = "\t", quote = FALSE, row.names = FALSE)
