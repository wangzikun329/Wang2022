# Read in the loadinginfo.csv, locate the DAM files and analyze sleep

library(data.table)
# ----------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

WORK_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
Batch <- "0722"
DATA_DIR <- "C:/Users/wangz/Desktop/Laboratory/SleepAnalysis/SleepLongevityData/DAMrawdata"
setwd(WORK_DIR)

loadinginfo <- read.csv(paste("loadinginfo_",Batch,".csv",sep = ""))
N_DAYS <- 5
monitorlist <- unique(loadinginfo$monitor)
genotypelist <- unique(loadinginfo$genotype)
# ----------------------------------------------------------------------------
# link monitor files
library(damr)
loadinginfo_linked <- link_dam_metadata(loadinginfo, result_dir = DATA_DIR)

#load in data, apply sleepr functions when loading
library(sleepr)
dt <- load_dam(loadinginfo_linked,FUN = sleepr::sleep_dam_annotation)
summary(dt)

# detecting anomalies
# plotting activity data based on monitors, save individual pictures to pdf files.
library(ggetho)
for(i in monitorlist)
{
  pdf(paste0(Batch,'_activity_by_monitor',i,'.pdf'))
  activity_by_monitor <-ggetho(dt[xmv(monitor) == i ], aes(z=activity)) +
    stat_bar_tile_etho() +
    scale_x_days()
  print(activity_by_monitor)
  dev.off()
}

# plotting sleep data based on monitors, save individual pictures.
for(i in monitorlist)
{
  pdf(paste0(Batch,'_sleep_by_monitor',i,'.pdf'))
  sleep_by_monitor <- ggetho(dt[xmv(monitor) == i], aes(z=asleep)) +
    stat_bar_tile_etho()+
    scale_x_days()
  print(sleep_by_monitor)
  dev.off()
}
# remove dead animals, generate a new behavior table called dt_curated
dt_curated_1 <- curate_dead_animals(dt)
summary(dt)
summary(dt_curated_1)
# find which ID is being removed in dt_curated
setdiff(dt[,id,meta=T],
        dt_curated_1[,id,meta=T])
curated_1_list<-data.table(setdiff(dt[,id,meta=T],
                                   dt_curated_1[,id,meta=T]))
fwrite(curated_1_list, paste("curated_1_list_",Batch,".csv",sep=""))
# plot sleep data before curation
pdf(paste0(Batch,'_sleep_before_deadcheck','.pdf'))
sleep_before_deadcheck<-ggetho(dt, aes(z=asleep)) +
  stat_ld_annotations(ypos = "top")+
  stat_tile_etho()
print(sleep_before_deadcheck)
dev.off()

# plot the sleep data after curation
pdf(paste0(Batch,'_sleep_after_deadcheck','.pdf'))
sleep_after_deadcheck<-ggetho(dt_curated_1, aes(z=asleep)) +
  stat_ld_annotations(ypos = "top")+
  stat_tile_etho()
print(sleep_after_deadcheck)
dev.off()
# further curation to remove animals dying too early
# we make a summary table of all lifespan for each animals
lifespan_dt <- dt_curated_1[, .(lifespan = max(t)), by=id]
# we filter this table for lifespan>4 and we keep the id
valid_ids <- lifespan_dt[lifespan > days(4), id]
invalid_ids <- lifespan_dt[lifespan <= days(4), id]
# we apply this filter
dt_curated_2 <- dt_curated_1[id %in% valid_ids]
summary(dt_curated_2)
# find which ID is being further removed
setdiff(dt_curated_1[,id,meta=T],
        dt_curated_2[,id,meta=T])
curated_2_list<-data.table(setdiff(dt_curated_1[,id,meta=T],
                                   dt_curated_2[,id,meta=T]))
fwrite(curated_2_list, paste("curated_2_list_",Batch,".csv"))

# mannually removing animal from analysis by modifying status column
# trimming data
# keeping data of only certain length
dt_curated_final <- dt_curated_2[t %between% c(days(0), days(N_DAYS))]

# population plots
pdf(paste0(Batch,'_population_sleep_1','.pdf'))
print(
  ggetho(dt_curated_final, aes(y=asleep, colour=genotype)) +
    stat_pop_etho() +
    stat_ld_annotations() +
    facet_grid(genotype ~ .) +
    scale_y_continuous(name= "Fraction of time sleeping",labels = scales::percent)
)
dev.off()
# wrap for a day, population level
pdf(paste0(Batch,'_population_sleep_wrap','.pdf'))
print(
  ggetho(dt_curated_final, aes(y=asleep, colour=genotype), time_wrap = hours(24)) +
    stat_pop_etho() +
    stat_ld_annotations() +
    facet_grid(genotype ~ .) +
    scale_y_continuous(name= "Fraction of time sleeping",labels = scales::percent)
)
dev.off()

# sleep profile
order <- c("WT"="black","inc1"="red","wakeD2"="green")
pdf(paste0(Batch,'_population_sleep_profile','.pdf'))
print(
  ggetho(dt_curated_final, aes(y=asleep*100, colour=genotype),time_wrap = hours(24)) +
    stat_pop_etho() +
    scale_y_continuous(name="Sleep(%)",limits=c(0,100),expand =c(0,0)) +
    scale_colour_manual(values = order, breaks = genotypelist,
                        labels = genotypelist,
                        aesthetics = c("colour","fill"))+
    theme_classic()+
    theme(
      aspect.ratio=0.4,
      axis.title.x = element_text(color = "black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "black", size = 14, face = "bold"),
      axis.text = element_text(colour = "black", size =14, face = "bold"),
      axis.ticks.length.y = unit(.25, "cm"),
      axis.ticks.length.x = unit(.25, "cm"),
      axis.text.x = element_text(margin = margin(t = 1, unit = "cm")),
      panel.grid.major.y = element_line(colour="gray90", size = (0.5)),
      legend.title = element_blank(),
      legend.text = element_text(colour="black", size=14, face="bold"),
      legend.position="top"
    )
)
dev.off()
# data summary
# day night sleep calculation using modulo operation
dt_curated_final[, phase := ifelse(t %% hours(24) < hours(12), "L", "D")]
summary_dt_final <- rejoin(dt_curated_final[,
                          .(
                            # this is where the computation happens
                            sleep_fraction_all = mean(asleep),
                            sleep_time_all = 1440*mean(asleep),
                            sleep_fraction_l = mean(asleep[phase == "L"]),
                            sleep_time_l = 720*mean(asleep[phase == "L"]),
                            sleep_fraction_d = mean(asleep[phase == "D"]),
                            sleep_time_d = 720*mean(asleep[phase == "D"])
                          ),
                          by=id])

#every 4hr bin  sleep calculations
dt_binned_sleep_4hr <- bin_apply_all(dt_curated_final, 
                                     asleep,
                                     x = t,
                                     x_bin_length = hours(4),
                                     wrap_x_by = days(1),
                                     FUN = mean
)

dt_sleep_ZT0_4 <- dt_binned_sleep_4hr[t == "0"]

summary_dt_sleep_ZT0_4 <- rejoin(dt_sleep_ZT0_4[,
                        .(
                          sleep_fraction_ZT0_4 = asleep,
                          sleep_time_ZT0_4 = 240*asleep
                        ),
                        by=id])

# preserve the id column, and the two columns giving ZT0-4 sleep parameters,
# combine this with the whole data table and write a file
summary_dt_sleep_ZT0_4[,c(1,13:14)]
summary_dt_final<-merge(summary_dt_final,summary_dt_sleep_ZT0_4[,c(1,13:14)],by="id")
summary_dt_final[,(c(2) ):= NULL]
fwrite(summary_dt_final,paste("summary_",Batch,".csv",sep=""))

# sleep architecture: bout analysis
bout_dt <- bout_analysis(asleep, dt_curated_final)
bout_dt <- bout_dt[asleep == TRUE, -"asleep"]
# bout length vs. time of the day
pdf(paste0(Batch, '_population_sleep_bout_wrap', '.pdf'))
print(
  ggetho(bout_dt, aes(y = duration / 60, colour=genotype), time_wrap = hours(24)) +
    stat_pop_etho() +
    facet_grid(genotype ~ .) +
    scale_y_continuous(name = "Bout length (min)")
)
dev.off()

# bout profile
pdf(paste0(Batch,'_population_sleep_bout_profile','.pdf'))
print(
  ggetho(bout_dt, aes(y = duration / 60, colour=genotype), time_wrap = hours(24)) +
    stat_pop_etho() +
    scale_y_continuous(name="Bout length (min)",expand =c(0,0)) +
    scale_colour_manual(values = order, breaks = genotypelist,
                        labels = genotypelist,
                        aesthetics = c("colour","fill"))+
    theme_classic()+
    theme(
      aspect.ratio=0.4,
      axis.title.x = element_text(color = "black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "black", size = 14, face = "bold"),
      axis.text = element_text(colour = "black", size =14, face = "bold"),
      axis.ticks.length.y = unit(.25, "cm"),
      axis.ticks.length.x = unit(.25, "cm"),
      axis.text.x = element_text(margin = margin(t = 1, unit = "cm")),
      panel.grid.major.y = element_line(colour="gray90", size = (0.5)),
      legend.title = element_blank(),
      legend.text = element_text(colour="black", size=14, face="bold"),
      legend.position="top"
    )
)
dev.off()

# latency to sleep_Day1
bout_dt_first_day <- bout_dt[t %between%  c(hours(0), hours(12))]
bout_summary <- rejoin(bout_dt_first_day[, .(
  Day1_latency = t[1],
  # the first bout is at t[1]
  Day1_first_bout_length = duration[1],
  Day1_latency_to_longest_bout = t[which.max(duration)]
),
by = id])
summary_dt_final<-merge(summary_dt_final,bout_summary[,c(1,13:15)],by="id")
# latency to sleep_Day2
bout_dt_second_day <- bout_dt[t %between%  c(days(1), days(1)+hours(12))]
bout_dt_second_day[, t:= t - days(1)]
bout_summary <- rejoin(bout_dt_second_day[, .(
  Day2_latency = t[1],
  # the first bout is at t[1]
  Day2_first_bout_length = duration[1],
  Day2_latency_to_longest_bout = t[which.max(duration)]
),
by = id])
summary_dt_final<-merge(summary_dt_final,bout_summary[,c(1,13:15)],by="id")
fwrite(summary_dt_final,paste("summary_with_latency_",Batch,".csv",sep=""))

# Call external code to calculate sleep parameters, analyze bout distribution
# and generate raster plot of sleep profile
CODE_PATH <- WORK_DIR
source(paste(CODE_PATH,"/Sleep_Parameters.R",sep = ""))
sleep_parameters(Batch)

source(paste(CODE_PATH,"/Sleep_Bout_Distribution.R",sep = ""))
bout_distribution(dt_curated_final, genotypelist, Batch)

source(paste(CODE_PATH,"/Sleep_Raster.R",sep = ""))
raster_plot(dt_curated_final,genotypelist,Batch)



