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
#Save all the bouts
#saveRDS(All_Bouts, file = "All_Bouts_DREADD_Batch1.rds")

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch1/Output")
l1 <- readRDS(file = "All_Bouts_taCasp3_Batch1.rds")
l2 <- readRDS(file = "All_Bouts_taCasp3_Batch2.rds")

#Merge the lists together
All_Bouts <- append(l1, l2)
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
Comparisons <- c("NotaCasp3 vs CoffFon-taCasp3",
                 "NotaCasp3 vs ConFon-taCasp3",
                 "CoffFon-taCasp3 vs ConFon-taCasp3")
ks_summary <- data.frame(matrix(nrow = 3*length(Parameters), ncol = 4))
colnames(ks_summary) <- c("Comparison", "Parameter", "D_statistic", "pvalue")
ks_summary$Comparison <- rep(Comparisons, times = length(Parameters))
ks_summary$Parameter <- rep(Parameters, each = length(Comparisons))

for (Row_n in 1:nrow(ks_summary)) {
  
  if(ks_summary[Row_n,]$Comparison == "NotaCasp3 vs CoffFon-taCasp3") {
    ks_result <- ks.test(x = df[df$Group == "NotaCasp3", ks_summary[Row_n,]$Parameter],
                         y = df[df$Group == "CoffFontaCasp3", ks_summary[Row_n,]$Parameter])
    
    ks_summary[Row_n,]$D_statistic <- ks_result$statistic
    ks_summary[Row_n,]$pvalue <- ks_result$p.value
  }
  
  if(ks_summary[Row_n,]$Comparison == "NotaCasp3 vs ConFon-taCasp3") {
    ks_result <- ks.test(x = df[df$Group == "NotaCasp3", ks_summary[Row_n,]$Parameter],
                         y = df[df$Group == "ConFontaCasp3", ks_summary[Row_n,]$Parameter])
    
    ks_summary[Row_n,]$D_statistic <- ks_result$statistic
    ks_summary[Row_n,]$pvalue <- ks_result$p.value
  }
  
  if(ks_summary[Row_n,]$Comparison == "CoffFon-taCasp3 vs ConFon-taCasp3") {
    ks_result <- ks.test(x = df[df$Group == "CoffFontaCasp3", ks_summary[Row_n,]$Parameter],
                         y = df[df$Group == "ConFontaCasp3", ks_summary[Row_n,]$Parameter])
    
    ks_summary[Row_n,]$D_statistic <- ks_result$statistic
    ks_summary[Row_n,]$pvalue <- ks_result$p.value
  }
  
}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch1/Output")
write.csv(ks_summary, "Bilateral-taCasp3_Bouts_Kolmogorov-Smirnov.csv")

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

Bouts_Freq <- Bouts_Freq %>%
  mutate(Sex = case_when(
    MouseID == "E701" ~ "M",
    MouseID == "E768" ~ "M",
    MouseID == "D696" ~ "F",
    MouseID == "E774" ~ "M",
    MouseID == "E772" ~ "F",
    MouseID == "E765" ~ "F",
    MouseID == "E767" ~ "M",
    MouseID == "E801" ~ "M",
    MouseID == "E770" ~ "F",
    MouseID == "F191" ~ "M",
    MouseID == "F305" ~ "M",
    MouseID == "F272" ~ "F",
    MouseID == "F285" ~ "F",
    MouseID == "F286" ~ "F",
    MouseID == "F190" ~ "M",
    MouseID == "F278" ~ "M",
    MouseID == "F482" ~ "M",
    MouseID == "F274" ~ "F",
    MouseID == "F276" ~ "F",
    MouseID == "F283" ~ "M",
    MouseID == "F290" ~ "M",
    MouseID == "F273" ~ "F",
    MouseID == "F275" ~ "F",
    MouseID == "F284" ~ "F"
  ))

Toplot <- Bouts_Freq
Toplot$Group <- factor(Toplot$Group, levels = c("NotaCasp3", "CoffFontaCasp3", "ConFontaCasp3"))

Toplot |>
  tidyplot(x = Group, y = Bouts_Freq, color = Group) |>
  adjust_colors(c("NotaCasp3" = "#B2B2B2", "CoffFontaCasp3" = "#FDB3B3", "ConFontaCasp3" = "#B2DCB4")) |> 
  add_boxplot() |> 
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  adjust_y_axis_title("Bouts Frequency") |>
  adjust_x_axis_title("")

#Perform one-way anova

#Summary statictic
Toplot %>%
  group_by(Group) %>%
  get_summary_stats(Bouts_Freq, type = "mean_sd")

#Test for outliers
Toplot %>%
  group_by(Group) %>%
  identify_outliers(Bouts_Freq)
#No extreme outliers were found

#Test for normality
model  <- lm(Bouts_Freq ~ Group, data = Toplot) # Build the linear model
ggqqplot(residuals(model)) # Create a QQ plot of residuals
shapiro_test(residuals(model)) # Compute saphiro-wilk
#Normality can be assumed

#Test for normality per group
Toplot %>%
  group_by(Group) %>%
  shapiro_test(Bouts_Freq)
#Note: The pValue is not significant so normality can be assumed

#Test homogeneity of variance
Toplot %>% levene_test(Bouts_Freq ~ Group)
#pValue is not significant so homogeneity of variance can be assumed

#ANOVA
res.aov <- Toplot %>% anova_test(Bouts_Freq ~ Group)
res.aov
#pValue was significant so I will proceed to Tukey post-hoc test

pwc <- Toplot %>% 
  tukey_hsd(Bouts_Freq ~ Group)
pwc

write.csv(pwc, "BoutsFreq_taCasp3_Bilateral_onewayAOV_TukeyHSD.csv")
