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

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/DREADD/Batch1/Output")
l1 <- readRDS(file = "All_Bouts_DREADD_Batch1.rds")
l2 <- readRDS(file = "All_Bouts_DREADD_Batch2.rds")

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

temp <- All_Bouts[["CoffFon-DREADD_I0252"]][[1]]
temp[temp$Frame > 0 & temp$Speed == 0, ]$Frame[1]



df <- data.frame(BoutID = BoutID,
           mouseID = mouseID,
           Group = Group,
           peakspeed = peakspeed,
           boutlenght = boutlenght,
           peakacceleration = peakacceleration,
           peakdeceleration = peakdeceleration 
           )


Parameters <- colnames(df[,4:7])
Comparisons <- c("fDIO-mCherry vs CoffFon-DREADD",
                 "fDIO-mCherry vs ConFon-DREADD",
                 "CoffFon-DREADD vs ConFon-DREADD")
ks_summary <- data.frame(matrix(nrow = 3*length(Parameters), ncol = 4))
colnames(ks_summary) <- c("Comparison", "Parameter", "D_statistic", "pvalue")
ks_summary$Comparison <- rep(Comparisons, times = length(Parameters))
ks_summary$Parameter <- rep(Parameters, each = length(Comparisons))

for (Row_n in 1:nrow(ks_summary)) {
  
  if(ks_summary[Row_n,]$Comparison == "fDIO-mCherry vs CoffFon-DREADD") {
    ks_result <- ks.test(x = df[df$Group == "fDIO-mCherry", ks_summary[Row_n,]$Parameter],
                         y = df[df$Group == "CoffFon-DREADD", ks_summary[Row_n,]$Parameter])
    
    ks_summary[Row_n,]$D_statistic <- ks_result$statistic
    ks_summary[Row_n,]$pvalue <- ks_result$p.value
  }
  
  if(ks_summary[Row_n,]$Comparison == "fDIO-mCherry vs ConFon-DREADD") {
    ks_result <- ks.test(x = df[df$Group == "fDIO-mCherry", ks_summary[Row_n,]$Parameter],
                         y = df[df$Group == "ConFon-DREADD", ks_summary[Row_n,]$Parameter])
    
    ks_summary[Row_n,]$D_statistic <- ks_result$statistic
    ks_summary[Row_n,]$pvalue <- ks_result$p.value
  }
  
  if(ks_summary[Row_n,]$Comparison == "CoffFon-DREADD vs ConFon-DREADD") {
    ks_result <- ks.test(x = df[df$Group == "CoffFon-DREADD", ks_summary[Row_n,]$Parameter],
                         y = df[df$Group == "ConFon-DREADD", ks_summary[Row_n,]$Parameter])
    
    ks_summary[Row_n,]$D_statistic <- ks_result$statistic
    ks_summary[Row_n,]$pvalue <- ks_result$p.value
  }
  
}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/DREADD/Batch1/Output")
write.csv(ks_summary, "DREADD_Bouts_Kolmogorov-Smirnov.csv")

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
    MouseID == "I0251" ~ "M",
    MouseID == "I0258" ~ "M",
    MouseID == "I0540" ~ "M",
    MouseID == "I0603" ~ "M",
    MouseID == "I0617" ~ "F",
    MouseID == "I0413" ~ "F",
    MouseID == "I0252" ~ "M",
    MouseID == "I0260" ~ "M",
    MouseID == "I0544" ~ "M",
    MouseID == "I0604" ~ "M",
    MouseID == "I0418" ~ "M",
    MouseID == "I0804" ~ "F",
    MouseID == "I0253" ~ "M",
    MouseID == "I0263" ~ "M",
    MouseID == "I0606" ~ "M",
    MouseID == "I0258" ~ "M",
    MouseID == "I0540" ~ "M",
    MouseID == "I0603" ~ "M",
    MouseID == "I0617" ~ "F",
    MouseID == "I0413" ~ "F",
    MouseID == "I0408" ~ "F",
    MouseID == "I0619" ~ "F",
    MouseID == "I0825" ~ "F",
    MouseID == "I0826" ~ "F",
    MouseID == "I0827" ~ "F",
    MouseID == "I0917" ~ "F",
    MouseID == "I0918" ~ "F",
    MouseID == "J0010" ~ "F",
    MouseID == "J0012" ~ "F",
    MouseID == "J0013" ~ "F",
    MouseID == "J0015" ~ "M"
  ))

Toplot <- Bouts_Freq
Toplot$Group <- factor(Toplot$Group, levels = c("fDIO-mCherry", "CoffFon-DREADD", "ConFon-DREADD"))

Toplot |>
  tidyplot(x = Group, y = Bouts_Freq, color = Group) |>
  adjust_colors(c("fDIO-mCherry" = "#B2B2B2", "CoffFon-DREADD" = "#FDB3B3", "ConFon-DREADD" = "#B2DCB4")) |> 
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

write.csv(pwc, "BoutsFreq_DREADD_onewayAOV_TukeyHSD.csv")
