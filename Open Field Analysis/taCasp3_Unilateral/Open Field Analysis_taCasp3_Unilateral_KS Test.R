library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(DescTools)
library(magick)
library(imager)
library(lattice)
library(viridis)
#library(paletteer)
library(gplots)
library(sommer)
library(tidyplots)
library(ggpubr)
library(rstatix)

#------------------------------------------------------------------------

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1/Output")
l1 <- readRDS(file = "All_Bouts_taCasp3_Batch1.rds")
l2 <- readRDS(file = "All_Bouts_taCasp3_Batch2.rds")
l3 <- readRDS(file = "All_Bouts_taCasp3_Batch3.rds")

#Merge the lists together
All_Bouts <- append(l1, l2)
All_Bouts <- append(All_Bouts, l3)

#------------------------------------------------------------------------
#Extract the peak speed, acceleration and deceleration for every bout
BoutID <- c()
mouseID <- c()
Group <- c()
peakspeed <- c()
boutlenght <- c()
peakacceleration <- c()
peakdeceleration <- c()

counter <- 1
for (Mouse_n in names(All_Bouts)) {
  
  for (Bout_n in 1:length(All_Bouts[[Mouse_n]])) {
    
    BoutID[counter] <- counter
    mouseID[counter] <- sapply(strsplit(Mouse_n, split = "_"), "[", 2)
    Group[counter] <- sapply(strsplit(Mouse_n, split = "_"), "[", 1)
    peakspeed[counter] <- max(All_Bouts[[Mouse_n]][[Bout_n]]$Speed)
    peakacceleration[counter] <- max(All_Bouts[[Mouse_n]][[Bout_n]]$Acceleration)
    peakdeceleration[counter] <- min(All_Bouts[[Mouse_n]][[Bout_n]]$Acceleration)
    
    temp <- All_Bouts[[Mouse_n]][[Bout_n]]
    boutlenght[counter] <- temp[temp$Frame > 0 & temp$Speed == 0, ]$Frame[1]
    
    counter <- counter+1
  }
}


df <- data.frame(BoutID = BoutID,
           mouseID = mouseID,
           Group = Group,
           peakspeed = peakspeed,
           boutlenght = boutlenght,
           peakacceleration = peakacceleration,
           peakdeceleration = peakdeceleration 
           )


Parameters <- colnames(df[,4:7])
Comparisons <- c("NotaCasp3 vs CoffFon-taCasp3")
ks_summary <- data.frame(matrix(nrow = length(Parameters), ncol = 4))
colnames(ks_summary) <- c("Comparison", "Parameter", "D_statistic", "pvalue")
ks_summary$Comparison <- rep(Comparisons, times = length(Parameters))
ks_summary$Parameter <- rep(Parameters, each = length(Comparisons))

for (Row_n in 1:nrow(ks_summary)) {
  
  ks_result <- ks.test(x = df[df$Group == "NotaCasp3", ks_summary[Row_n,]$Parameter],
                       y = df[df$Group == "CoffFontaCasp3", ks_summary[Row_n,]$Parameter])
  ks_summary[Row_n,]$D_statistic <- ks_result$statistic
  ks_summary[Row_n,]$pvalue <- ks_result$p.value
 
} 
  

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1/Output")
write.csv(ks_summary, "Unilateral-taCasp3_Bouts_Kolmogorov-Smirnov.csv")

#------------------------------------------------------------------------
#Plot bouts frequency
nminutes <- 15 #the velocity was recorded for 15 minutes

Bouts_Freq <- as.data.frame(matrix(nrow = length(All_Bouts), ncol = 4))
colnames(Bouts_Freq) <- c("MouseID", "Group", "Sex", "Bouts_Freq")
Bouts_Freq$MouseID <- sapply(strsplit(names(All_Bouts), "_"), "[", 2)
Bouts_Freq$Group <- sapply(strsplit(names(All_Bouts), "_"), "[", 1)
rownames(Bouts_Freq) <- names(All_Bouts)

for (MouseID in names(All_Bouts)) {
  Bouts_Freq[MouseID, "Bouts_Freq"] <- length(All_Bouts[[MouseID]]) / nminutes
}         
Bouts_Freq
Bouts_Freq <- Bouts_Freq %>%
  mutate(Sex = case_when(
    MouseID == "F574" ~ "M",
    MouseID == "F707" ~ "M",
    MouseID == "F658" ~ "M",
    MouseID == "F708" ~ "M",
    MouseID == "F574" ~ "M",
    MouseID == "F706" ~ "F",
    MouseID == "G0262" ~ "M",
    MouseID == "G0215" ~ "F",
    MouseID == "G0217" ~ "F",
    MouseID == "G0350" ~ "F",
    MouseID == "G0264" ~ "M",
    MouseID == "G0218" ~ "F",
    MouseID == "G0219" ~ "F",
    MouseID == "G0351" ~ "F",
    MouseID == "I0135" ~ "M",
    MouseID == "I0137" ~ "M",
    MouseID == "I0214" ~ "M",
    MouseID == "I0029" ~ "F",
    MouseID == "I0031" ~ "F",
    MouseID == "I0250" ~ "F",
    MouseID == "I0138" ~ "M",
    MouseID == "I0212" ~ "M",
    MouseID == "I0213" ~ "M",
    MouseID == "I0032" ~ "F",
    MouseID == "I0248" ~ "F",
    MouseID == "I0256" ~ "F"
  ))

Toplot <- Bouts_Freq
Toplot$Group <- factor(Toplot$Group, levels = c("NotaCasp3", "CoffFontaCasp3"))

Toplot |>
  tidyplot(x = Group, y = Bouts_Freq, color = Group) |>
  adjust_colors(c("NotaCasp3" = "#B2B2B2", "CoffFontaCasp3" = "#FDB3B3")) |> 
  add_boxplot() |> 
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  adjust_y_axis_title("Bouts Frequency") |>
  add_test_asterisks(hide_info = TRUE, method = "t.test") |> #Perform T-test
  adjust_x_axis_title("")

