
bout_distribution <- function(dt_curated_final, genotypelist, Batch=""){
  
  # sleep architecture: bout analysis
  bout_dt <- bout_analysis(asleep, dt_curated_final)
  # add day/night information 
  bout_dt <- bout_dt[, phase := ifelse(t %% hours(24) < hours(12), "L", "D")]
  # change duratioin to minute
  bout_dt_min <- bout_dt[,duration:=duration/60]
  # change t to minute 
  bout_dt_min <- bout_dt_min[,t:=t/60]
  # add lifespan column
  bout_dt_min <- bout_dt_min[, t_days:= t/1440]
  # remove asleep column
  bout_dt_min<-bout_dt_min[asleep == TRUE, -"asleep"]
  bout_dt_min<-bout_dt_min[duration>=5]
  # bout number and bout length
  summary_bout_dt <- rejoin(bout_dt_min[,
                                        .(
                                          # this is where the computation happens
                                          bout_number = length(duration)/max(t_days),
                                          bout_length = mean(duration),
                                          bout_number_l = length(duration[phase == "L"])/max(t_days),
                                          bout_length_l = mean(duration[phase == "L"]),
                                          bout_number_d = length(duration[phase == "D"])/max(t_days),
                                          bout_length_d = mean(duration[phase == "D"])
                                        ),
                                        by=id])
  summary_bout_dt[,(c(2) ):= NULL]
  fwrite(summary_bout_dt,paste("summary_sleep_bouts_",Batch,".csv",sep=""))
  
  # calculate mean, std, N of sleep bouts for each genotype
  bout_statistics <- data.frame(matrix(data=NA,nrow=length(genotypelist),ncol=18),
                                row.names = genotypelist)
  names(bout_statistics) <- paste(rep(names(summary_bout_dt)[12:17],each=3),c("mean","sd","N"),sep="_")
  for(i in genotypelist){
    bout_data <- summary_bout_dt[summary_bout_dt[,genotype == i]]
    bout_statistics[i,"bout_number_mean"] = mean(bout_data$bout_number)
    bout_statistics[i,"bout_number_sd"] = sd(bout_data$bout_number,na.rm = TRUE)
    bout_statistics[i,"bout_number_N"] = length(bout_data$bout_number)
    bout_statistics[i,"bout_number_l_mean"] = mean(bout_data$bout_number_l)
    bout_statistics[i,"bout_number_l_sd"] = sd(bout_data$bout_number_l,na.rm = TRUE)
    bout_statistics[i,"bout_number_l_N"] = length(bout_data$bout_number_l)
    bout_statistics[i,"bout_number_d_mean"] = mean(bout_data$bout_number_d)
    bout_statistics[i,"bout_number_d_sd"] = sd(bout_data$bout_number_d,na.rm = TRUE)
    bout_statistics[i,"bout_number_d_N"] = length(bout_data$bout_number_d)
    bout_statistics[i,"bout_length_mean"] = mean(bout_data$bout_length)
    bout_statistics[i,"bout_length_sd"] = sd(bout_data$bout_length,na.rm = TRUE)
    bout_statistics[i,"bout_length_N"] = length(bout_data$bout_length)
    bout_statistics[i,"bout_length_l_mean"] = mean(bout_data$bout_length_l)
    bout_statistics[i,"bout_length_l_sd"] = sd(bout_data$bout_length_l,na.rm = TRUE)
    bout_statistics[i,"bout_length_l_N"] = length(bout_data$bout_length_l)
    bout_statistics[i,"bout_length_d_mean"] = mean(bout_data$bout_length_d)
    bout_statistics[i,"bout_length_d_sd"] = sd(bout_data$bout_length_d,na.rm = TRUE)
    bout_statistics[i,"bout_length_d_N"] = length(bout_data$bout_length_d)
  }
  bout_statistics["Genotype"] <- paste(row.names(bout_statistics),Batch,sep = "_")
  bout_statistics <- bout_statistics[,c(19,1:18)]
  fwrite(bout_statistics,paste(Batch,"_sleep_bout_statistics.csv",sep=""))

  # bout distrubution--------------------------------------------------------
  bout_dt_min_L<-bout_dt_min[phase=="L"]
  bout_dt_min_D<-bout_dt_min[phase=="D"]
  
  # daytime sleep bout histogram/density plots, 
  # pay attention to argument origin and right, opposite to hist function, right=TRUE means right-open, [a,b)
  order <- c("WT"="black","inc1"="red","wakeD2"="green")
  # b<-bout_dt_min_L[xmv(genotype)==genotypelist[1]|xmv(genotype)==genotypelist[2]||xmv(genotype)==genotypelist[3]]
  # b<-b[,genotype:=xmv(genotype)]
  bout_dt_min_L <- bout_dt_min_L[,genotype:=xmv(genotype)]
  pdf(paste0(Batch,'_DayTime_SleepBouts','.pdf'))
  print(
    ggplot(bout_dt_min_L, aes(x=duration, colour=genotype, fill=genotype)) +
      geom_histogram(aes(y=..density..), alpha=0.5, binwidth=10, origin=5, right=TRUE,
                     position="identity")+
      scale_colour_manual(values = order, breaks = genotypelist,
                          labels = genotypelist,
                          aesthetics = c("colour","fill"))+
      facet_grid(genotype ~ .)+
      scale_x_continuous(name="Duration of Daytime Sleep Bout (min)", limit=c(0,720)) +
      scale_y_continuous(limit=c(0,0.065)) +
      theme_classic()+
      theme(
        aspect.ratio=1,
        axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        axis.text = element_text(colour = "black", size =14, face = "bold"),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(.25, "cm"),
        panel.grid.major.y = element_line(colour="gray90", size = (0.5)),
        axis.text.x = element_text(margin = margin(t = 0.25, unit = "cm")),
        axis.text.y = element_text(margin = margin(t = 0.25, unit = "cm")),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=14, face="bold"),
        strip.text.y = element_blank(),
        legend.position="top"
      )
  )
  dev.off()
  
  
  bout_dt_min_D <- bout_dt_min_D[,genotype:=xmv(genotype)]
  pdf(paste0(Batch,'_NightTime_SleepBouts','.pdf'))
  print(
    ggplot(bout_dt_min_D, aes(x=duration, colour=genotype, fill=genotype)) +
      geom_histogram(aes(y=..density..), alpha=0.5, binwidth=10, origin=5, right=TRUE,
                     position="identity")+
      scale_colour_manual(values = order, breaks = genotypelist,
                          labels = genotypelist,
                          aesthetics = c("colour","fill"))+
      facet_grid(genotype ~ .)+
      scale_x_continuous(name="Duration of Daytime Sleep Bout (min)", limit=c(0,720)) +
      scale_y_continuous(limit=c(0,0.065)) +
      theme_classic()+
      theme(
        aspect.ratio=1,
        axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        axis.text = element_text(colour = "black", size =14, face = "bold"),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(.25, "cm"),
        panel.grid.major.y = element_line(colour="gray90", size = (0.5)),
        axis.text.x = element_text(margin = margin(t = 0.25, unit = "cm")),
        axis.text.y = element_text(margin = margin(t = 0.25, unit = "cm")),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=14, face="bold"),
        strip.text.y = element_blank(),
        legend.position="top"
      )
  )
  dev.off()
  
  order <- c("WT"=1,"inc1"=2,"wakeD2"=3)  
  # take genotype, subsetting the bout_min_l or bout_min_d table, make frequency counts, write to a table
  bout_l_dist <- data.frame("bout_day"=c(1:1499))
  
  for (i in genotypelist)
  {
    a<-bout_dt_min_L[xmv(genotype)==i]
    amax<-max(a[,duration])
    factor<-factor(a[,duration],levels=1:amax)
    out <- as.data.frame(table(factor))
    out <- transform(out, cumFreq = cumsum(Freq), prob = prop.table(Freq))
    out <- transform(out, cumProp = cumsum(prob))
    out<-tibble::rownames_to_column(out)
    names(out) <- paste(names(out),order[i],i,Batch,sep = "_")
    out <- out[,2:6]
    names(out)[1] <- "bout_day"
    bout_l_dist <- merge(bout_l_dist, out, all = TRUE)
    #fwrite(out, paste(Batch,"_DayTimeDistribution_",i,".csv",sep=""))
  }
  bout_l_dist <- bout_l_dist[,sort(names(bout_l_dist))] 
  fwrite(bout_l_dist, paste(Batch,"_DayTimeDistribution.csv",sep=""))
  

  bout_d_dist <- data.frame("bout_night"=c(1:1499))
  for (i in genotypelist)
  {
    a<-bout_dt_min_D[xmv(genotype)==i]
    amax<-max(a[,duration])
    factor<-factor(a[,duration],levels=1:amax)
    out <- as.data.frame(table(factor))
    out <- transform(out, cumFreq = cumsum(Freq), prob = prop.table(Freq))
    out <- transform(out, cumProp = cumsum(prob))
    out<-tibble::rownames_to_column(out)
    names(out) <- paste(names(out),order[i],i,Batch,sep = "_")
    out <- out[,2:6]
    names(out)[1] <- "bout_night"
    bout_d_dist <- merge(bout_d_dist, out, all = TRUE)
    #fwrite(out, paste(Batch,"_NightTimeDistribution_",i,".csv",sep=""))
  }
  bout_d_dist <- bout_d_dist[,sort(names(bout_d_dist))] 
  fwrite(bout_d_dist, paste(Batch,"_NightTimeDistribution.csv",sep=""))
  
  #take genotype, subsetting the bout_min_l or bout_min_d table, write to a table of bout lists
  for (i in genotypelist)
  {
    a<-bout_dt_min_L[xmv(genotype)==i]
    boutlist<-as.data.table(a[,duration])
    names(boutlist)[1]<-paste0(i,"_DayTimeBouts")
    fwrite(boutlist, paste(Batch,"_DayTimeBouts_",i,".csv",sep=""))
  }

  for (i in genotypelist)
  {
    a<-bout_dt_min_D[xmv(genotype)==i]
    boutlist<-as.data.table(a[,duration])
    names(boutlist)[1]<-paste0(i,"_NightTimeBouts")
    fwrite(boutlist, paste(Batch,"_NightTimeBouts_",i,".csv",sep=""))
  }

}



