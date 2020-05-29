# write loading information
# For every experiment, manualy change <Batch>, <monitorlist>,
# <genotypelist>, and <loadinginfo>

library(data.table)
# ----------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

DATA_DIR <- "C:/Users/wangz/Desktop/Laboratory/SleepAnalysis/SleepLongevityData/DAMrawdata"
WORK_DIR <- paste("C:/Users/wangz/Desktop/Laboratory/SleepAnalysis/SleepLongevityData/", Batch, sep = "")
setwd(WORK_DIR)

Batch = "0722"
monitorlist <- c("M39","M51","M57")
genotypelist <- c("WT","wakeD2","inc1")
# ----------------------------------------------------------------------------
# write experiment design
library(data.table)
loadinginfo <- data.table(
  file = rep(c("Monitor39.txt", "Monitor51.txt", "Monitor57.txt"), each = 32),
  monitor = rep(monitorlist, each = 32 ),
  region_id = 1:32,
  status = "OK",
  start_datetime = "2016-07-23 11:00:00",
  stop_datetime = "2016-07-28 12:00:00",
  sex = "M",
  genotype = rep(genotypelist, each = 32),
  BDSC = rep(genotypelist, each = 32),
  treatment = rep(c("NA","NA","NA"),each = 32)
)

# write loadinginfo file
fwrite(loadinginfo, paste("loadinginfo_",Batch,".csv",sep = ""))
