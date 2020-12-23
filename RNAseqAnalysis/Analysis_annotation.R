
# the output file of this script could not be read by R
file.name <- paste0("DGE_TC_FINAL",".txt")

# install.packages("remotes")
# remotes::install_github("hangnoh/flybaseR")

# setup
# setwd("~/Documents/Laboratory/RNA_seq")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(org.Dm.eg.db)

# read in file with Flybase ID as row names
read.file <- read.table(file = file.name, header = TRUE)
# id.current <- read.csv(file = "FLYBASE_conversion.csv", header = TRUE, row.names = 1)

# retrieve annotation from database
# id.new <- id.current[row.names(read.file),]
# anno.DGE <- select(org.Dm.eg.db, keys = as.character(id.new$current_ID), keytype = "FLYBASE", columns=c("FLYBASECG","GENENAME","ENTREZID","ENSEMBL"))
anno.DGE <- select(org.Dm.eg.db, keys = row.names(read.file), keytype = "FLYBASE", columns=c("FLYBASECG","SYMBOL","GENENAME","ENTREZID"))
out.df <- cbind(anno.DGE, read.file)

write.table(out.df, file = paste(file.name,"anno.txt", sep="."), sep = "\t", quote = FALSE, row.names = FALSE)
