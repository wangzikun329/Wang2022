
setwd("~/Documents/Laboratory/RNA_seq")
library(NMF)
library(ggplot2)
library(circular)
library(tidyverse)

setwd("~/Documents/Laboratory/RNA_seq/")

# BiocManager::install("tidyverse")

wt2u.union <- read.table(file = "circadian_union_JTK_WT2U_result.txt", header = TRUE)
inc1.union <- read.table(file = "circadian_union_JTK_INC1_result.txt", header = TRUE)
waked2.union <- read.table(file = "circadian_union_JTK_WAKED2_result.txt", header = TRUE)
wt2uold.union <- read.table(file = "circadian_union_JTK_WT2U_old_result.txt", header = TRUE)

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
    annotate("text", x = rep(1, 7), y = c(0, 50, 100, 150, 200, 250, 300), label = c("0", "50", "100", "150", "200", "250", "300"), color="black",alpha=0.6, size=2, angle=0, fontface="bold", hjust=1) +
    ylim(-100,300) +
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
    geom_text(data=label_data, aes(x=as.factor(mids-1), y=5, label=counts, hjust=hjust), color="black", fontface="bold",alpha=1, size=3, angle= label_data$angle, inherit.aes = TRUE ) 
  
  if(!is.null(savename)){
    hist.plot <- hist.plot + ggtitle(savename) + theme(plot.title = element_text(hjust = 0.5))
    pdf(paste(savename, ".pdf", sep = ""))
    print(hist.plot)
    dev.off()
  }
  hist.plot
}

circular_hist(wt2u.union)
circular_hist(inc1.union)
circular_hist(waked2.union)
circular_hist(wt2uold.union)

circular_hist(wt2u.union, savename = "circadian_phase_hist_WT2U")
circular_hist(inc1.union, savename = "circadian_phase_hist_INC1")
circular_hist(waked2.union, savename = "circadian_phase_hist_WAKED2")
circular_hist(wt2uold.union, savename = "circadian_phase_hist_WT2U_old")


# read in circadian count data and phase data
wt2u.jtk <- read.table(file = "JTK_CYCLE_WT2U_out.txt", header = TRUE, row.names = "x")
inc1.jtk <- read.table(file = "JTK_CYCLE_INC1_out.txt", header = TRUE, row.names = "x")
waked2.jtk <- read.table(file = "JTK_CYCLE_WAKED2_out.txt", header = TRUE, row.names = "x")
wt2uold.jtk <- read.table(file = "JTK_CYCLE_WT2U_old_out.txt", header = TRUE, row.names = "x")

inc1.wt2u <- inc1.jtk[row.names(wt2u.union),]
waked2.wt2u <- waked2.jtk[row.names(wt2u.union),]
wt2uold.wt2u <- wt2uold.jtk[row.names(wt2u.union),]

circular_hist(wt2u.union)
circular_hist(inc1.wt2u)
circular_hist(waked2.wt2u)
circular_hist(wt2uold.wt2u)


wt2u.wt2uold <- wt2u.jtk[row.names(wt2uold.union),]
circular_hist(wt2uold.union)
circular_hist(wt2u.wt2uold)
