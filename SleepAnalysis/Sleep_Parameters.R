# Calculate total sleep, daytime sleep, nighttime sleep
# ----------------------------------------------------------------------------

sleep_parameters <- function(Batch){
  readin_summary_dt_final <- fread(paste("summary_",Batch,".csv",sep = ""))
  
  
  sleep_time_all<-dcast(readin_summary_dt_final,id~genotype,value.var=c("sleep_time_all"))
  fwrite(sleep_time_all, paste("sleep_time_all_",Batch,".csv",sep = ""))
  
  sleep_time_l<-dcast(readin_summary_dt_final,id~genotype,value.var=c("sleep_time_l"))
  fwrite(sleep_time_l, paste("sleep_time_l_",Batch,".csv",sep = ""))
  
  sleep_time_d<-dcast(readin_summary_dt_final,id~genotype,value.var=c("sleep_time_d"))
  fwrite(sleep_time_d, paste("sleep_time_d_",Batch,".csv",sep = ""))
  
  sleep_time_ZT0_4<-dcast(readin_summary_dt_final,id~genotype,value.var=c("sleep_time_ZT0_4"))
  fwrite(sleep_time_ZT0_4, paste("sleep_time_ZT0_4_",Batch,".csv",sep = ""))
}


