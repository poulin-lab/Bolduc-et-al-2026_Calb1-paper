library(readxl)
library(tidyverse)
library(data.table)
library(tidyplots)

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Manuscripts/taCasp3")

#---------------------------------------------------------------------------------------
#Bilateral

data <- readxl::read_excel("taCasp3_WeightLoss.xlsx", sheet = "Feuil1")

MouseID <- rep(colnames(data)[-1], each = length(data$`Days elapsed`))
Days_Elapsed <- rep(data$`Days elapsed`, times = length(colnames(data)[-1]))
Weight <- as.vector(as.matrix(data[,2:ncol(data)]))

df <- data.frame( MouseID = MouseID, Days_Elapsed = Days_Elapsed, Weight = Weight)
df$Group <- sapply(strsplit(df$MouseID, split = "_"), `[`, 1)
df <- drop_na(df)

df |>
  tidyplot(x = Days_Elapsed, y = Weight, color = Group, dodge_width = 0) |>
  add_mean_line(group = MouseID, alpha = 0.5, linewidth = 1) |>
  adjust_y_axis_title("% of Initial Weight") |>
  adjust_y_axis_title("Days Post-Injection") |>
  adjust_x_axis(rotate_labels = 90, breaks = seq(from = 0, to = 28, by = 4)) 

#---------------------------------------------------------------------------------------
#Unilateral
df <- read.csv("taCasp3_WeightLoss_Uni.csv")
df <- drop_na(df)

df |>
  tidyplot(x = Days, y = WeightLoss, color = Group, dodge_width = 0) |>
  add_mean_line(group = MouseID, alpha = 0.5, linewidth = 1) |>
  adjust_y_axis_title("% of Initial Weight") |>
  adjust_y_axis_title("Days Post-Injection") |>
  adjust_x_axis(rotate_labels = 90, breaks = seq(from = 0, to = 28, by = 4)) 

