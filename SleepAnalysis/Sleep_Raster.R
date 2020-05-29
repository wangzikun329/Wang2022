#sleep raster plots
library("ggplot2")

raster_plot <- function(dt_curated_final, genotypelist, Batch=""){
  dt_plot<-dt_curated_final[t %between% c(days(0), days(2))]
  order <- c("WT"="black","inc1"="red","wakeD2"="green")
  
  for(i in genotypelist){
    pdf(paste0(Batch,'_SleepRaster_',i,'.pdf'))
    print(
      ggplot(dt_plot[xmv(genotype)==i],aes(x=t,y=id,fill=asleep,height=0.70)) +
        geom_tile () +
        scale_y_discrete("Fly#") +
        scale_fill_manual(values = c("#FFFFFF",	"#FF2600")) +
        theme_minimal()+
        scale_x_hours()+
        theme_classic()+
        theme(
          aspect.ratio=0.4,
          axis.title.x = element_text(color = "black", size = 14, face = "bold"),
          axis.title.y = element_text(color = "black", size = 14, face = "bold"),
          axis.text.y = element_blank(),
          axis.text = element_text(color = "black", size = 14, face = "bold"),
          axis.ticks.length.y = unit(.10, "cm"),
          axis.ticks.length.x = unit(.25, "cm"),
          axis.text.x = element_text(margin = margin(t = 1, unit = "cm")),
          panel.grid.major.y = element_line(colour="gray90", size = (0.25)),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none")
    )
    dev.off()
  }
  
}
