
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(NMF)
library(ggplot2)
library(circular)
library(tidyverse)

# BiocManager::install("tidyverse")

wt2u.jtk <- read.table(file = "JTK_CYCLE_WT2U_out.txt", header = TRUE, row.names = "x")
inc1.jtk <- read.table(file = "JTK_CYCLE_INC1_out.txt", header = TRUE, row.names = "x")
waked2.jtk <- read.table(file = "JTK_CYCLE_WAKED2_out.txt", header = TRUE, row.names = "x")
wtold.jtk <- read.table(file = "JTK_CYCLE_WT2U_old_out.txt", header = TRUE, row.names = "x")

inc1.diff <- read.table(file = "circadian_setdiff_INC1_out.txt", header = TRUE)
waked2.diff <- read.table(file = "circadian_setdiff_WAKED2_out.txt", header = TRUE)
wtold.diff <- read.table(file = "circadian_setdiff_WT2U_old_out.txt", header = TRUE)

inc1.loss <- wt2u.jtk[row.names(inc1.diff)[1:161],]
waked2.loss <- wt2u.jtk[row.names(waked2.diff)[1:153],]
wtold.loss <- wt2u.jtk[row.names(wtold.diff)[1:111],]
mutant.loss <- wt2u.jtk[intersect(row.names(inc1.loss), row.names(waked2.loss)),]

circular_hist <- function(data.union, savename=NULL){
  lag.data <- data.union$LAG
  lag.data <- ifelse(lag.data >=24, lag.data-24, lag.data)
  hist.data <- hist(lag.data, breaks = seq(0,24,by=2), right = FALSE, plot = FALSE)
  hist.df <- data.frame(mids = hist.data$mids, counts = hist.data$counts, density = hist.data$density)
  
  label_data <- hist.df
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360/number_of_bar * (label_data$mids-1)/2
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  hist.plot <- ggplot(hist.df, aes(x=as.factor(mids-1), y =counts)) +
    geom_col(position = "stack", fill = "skyblue", colour = "black", width = 1, linetype = 1) +
    annotate("text", x = rep(1, 4), y = c(0, 10, 20, 30), label = c("0", "10", "20", "30"), color="black",alpha=0.6, size=2, angle=0, fontface="bold", hjust=1) +
    ylim(-10,35) +
    theme_minimal() +
    theme(
    #   axis.text = element_blank(),
       axis.title = element_blank(),
       axis.text.y = element_blank(),
       panel.grid.major.y = element_blank()
    #   panel.grid.major = element_line(),
    #   # plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
    ) +
    coord_polar(start = -0.25) +
    geom_text(data=label_data, aes(x=as.factor(mids-1), y=1, label=counts, hjust=hjust), color="black", fontface="bold",alpha=1, size=3, angle= label_data$angle, inherit.aes = TRUE )
    # geom_text(data=label_data, aes(x=as.factor(mids-1), y=counts+1, label=counts, hjust=hjust), color="black", fontface="bold",alpha=1, size=3, angle= label_data$angle, inherit.aes = TRUE )
  
  if(!is.null(savename)){
    hist.plot <- hist.plot + ggtitle(savename) + theme(plot.title = element_text(hjust = 0.5))
    pdf(paste(savename, ".pdf", sep = ""))
    print(hist.plot)
    dev.off()
  }
  hist.plot
}

# circular_hist(wt2u.union)
circular_hist(inc1.loss)
circular_hist(waked2.loss)
circular_hist(wtold.loss)
circular_hist(mutant.loss)

circular_hist(inc1.loss, savename = "circadian_loss_phase_hist_INC1")
circular_hist(waked2.loss, savename = "circadian_loss_phase_hist_WAKED2")
circular_hist(wtold.loss, savename = "circadian_loss_phase_hist_WTOLD")
circular_hist(mutant.loss, savename = "circadian_loss_phase_hist_mutant")

# write.table(inc1.loss, file = "circadian_loss_INC1_stats.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(waked2.loss, file = "circadian_loss_WAKED2_stats.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(wtold.loss, file = "circadian_loss_WTold_stats.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(mutant.loss, file = "circadian_loss_mutant_stats.txt", sep = "\t", quote = FALSE, row.names = TRUE)

wtold.gain <- wt2u.jtk[row.names(wtold.diff)[112:314],]
circular_hist(wtold.gain)
circular_hist(wtold.loss, savename = "circadian_gain_phase_hist_WTOLD")
