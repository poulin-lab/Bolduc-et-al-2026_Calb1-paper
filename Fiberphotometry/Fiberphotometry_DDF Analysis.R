library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(tidyplots)


#-------------------------------------------------------------------------------
################################################################################
##########################  PROCESSING OF DFF TRACES  ##########################
################################################################################
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/DFF/LowBF3Hz_IRLS300sec")
#Import the traces as a list
csv_files <- list.files(pattern = "\\.csv$")

Traces <- list()
for (n in 1:length(csv_files)) {
  Traces[[n]] <- read.csv(csv_files[n])
}

csv_files <- sapply(strsplit(csv_files, split =  "_"), "[[", 1)
names(Traces) <- csv_files
#-------------------------------------------------------------------------------
#Export the traces as graph
#Filter the traces to be included only in the first 25 min
Filtered_Traces <- Traces
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]] <- Traces[[n]][Traces[[n]]$Time < 1500 & Traces[[n]]$Time > 5,]
}

#Add AIN01 which is the motion corrected DFF (%)
for (mouseID in names(Filtered_Traces)) {
  Filtered_Traces[[mouseID]]$AIN01 <- Filtered_Traces[[mouseID]]$DFF * 100
}

#Export the DDF dLight representative traces 
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/DFF/LowBF3Hz_IRLS300sec")
for (n in 1:length(csv_files)) {
  
  plot <- ggplot(Filtered_Traces[[n]], aes(x=Time, y=AIN01)) +
    geom_line() +
    xlim(0,1500) + xlab("Time (sec.)") + scale_x_continuous(breaks = seq(0, 1500, by = 300)) +
    ylim(-15, 15) + 
    ylab("dLight deltaF/F (%)") + 
    theme_classic()
  
  ggsave(paste(csv_files[n], "_Trace.pdf", sep = ""),
         width = 8, 
         height = 4, 
         plot = plot)
}

#Export the DDF dLight representative traces for 100 sec.
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/DFF")
for (n in 1:length(csv_files)) {
  
  plot <- ggplot(Filtered_Traces[[n]], aes(x=Time, y=AIN01)) +
    geom_line() +
    xlim(0,1500) + 
    xlab("Time (sec.)") + 
    #scale_x_continuous(breaks = 20) +
    ylim(-15, 15) + 
    ylab("dLight deltaF/F (%)") + 
    theme_classic()
  
  ggsave(paste(csv_files[n], "_100sec_Trace.pdf", sep = ""),
         width = 8, 
         height = 4, 
         plot = plot)
}

#-------------------------------------------------------------------------------
#Define the time and bin width for generating the CI95% of the dLight standard deviation
max_time <- 1500
bin_width <- 100

#Bin the DeltaF/F values
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]]$Bin <- cut(Filtered_Traces[[n]]$Time, breaks = seq(from = 0, to = max_time, by = bin_width))
  Filtered_Traces[[n]]$Bin <- as.numeric(Filtered_Traces[[n]]$Bin)
}

#Generate a dataframe to store the standard error per bin for each mouse
DFF_SD <- data.frame(matrix(NA, nrow = max_time/bin_width, ncol = length(csv_files)))
rownames(DFF_SD) <- unique(rownames(Filtered_Traces[[1]]$Bin))
colnames(DFF_SD) <- csv_files

for (Bin_x in 1:nrow(DFF_SD)) {
  for (sample_y in 1:ncol(DFF_SD)) {
    DFF_SD[Bin_x, sample_y] <- sd(Filtered_Traces[[sample_y]][Filtered_Traces[[sample_y]]$Bin == Bin_x, ]$AIN01)
  }
}

#Stack the dataframe
DFF_SD <- stack(DFF_SD)

#Add a Time column representing the actual bin time values
DFF_SD$Time <- rep(seq(from = bin_width, to = max_time, by = bin_width), times = length(csv_files))


#Add a group column
DFF_SD <- DFF_SD %>%
  mutate(
    Group = case_when(
      ind == "J0782" ~ "fDIO-mCherry",
      ind == "J0859" ~ "fDIO-mCherry",
      ind == "J0777" ~ "fDIO-mCherry",
      ind == "K0992" ~ "fDIO-mCherry",
      ind == "K0871" ~ "fDIO-mCherry",
      ind == "J0871" ~ "fDIO-mCherry",
      
      ind == "J0784" ~ "CoffFon-DREADD",
      ind == "J0861" ~ "CoffFon-DREADD",
      ind == "K0861" ~ "CoffFon-DREADD",
      ind == "J0779" ~ "CoffFon-DREADD",
      ind == "K0331" ~ "CoffFon-DREADD",
      ind == "K0994" ~ "CoffFon-DREADD",
      
      ind == "K0023" ~ "ConFon-DREADD",
      ind == "K0029" ~ "ConFon-DREADD",
      ind == "K0863" ~ "ConFon-DREADD",
      ind == "K0865" ~ "ConFon-DREADD",
      ind == "K0990" ~ "ConFon-DREADD",
      
    )
  )

#Plot the graph of dLight DFF SD
DFF_SD |>
  tidyplot(x = Time, y = values, color = Group, dodge_width = 0) |>
  add_mean_line() |>
  add_ci95_ribbon() |>
  adjust_x_axis_title("Time (sec.)") |>
  adjust_x_axis(breaks = seq(from = 0, to = 1500, by = 300)) |>
  adjust_y_axis_title("dLight deltaF/F Standard Deviation")

#-------------------------------------------------------------------------------
################################################################################
##########################  PROCESSING OF dLight TRACES  #######################
################################################################################
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/dLight")

#Import the traces as a list
csv_files <- list.files(pattern = "\\.csv$")

Traces <- list()
for (n in 1:length(csv_files)) {
  Traces[[n]] <- read.csv(csv_files[n])
}

names(Traces) <- csv_files

#-------------------------------------------------------------------------------
#Export the traces as graph
#Filter the traces to be included only in the first 25 min
Filtered_Traces <- Traces
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]] <- Traces[[n]][Traces[[n]]$Time < 1500 & Traces[[n]]$Time > 5,]
}

#Normalise the curve value on the average value of the first 100 seconds
for (mouseID in names(Filtered_Traces)) {
  Filtered_Traces[[mouseID]]$FoldChange <- Filtered_Traces[[mouseID]]$AIN01 /mean(Filtered_Traces[[mouseID]][Filtered_Traces[[mouseID]]$Time < 100,]$AIN01)

}

#Export the DDF dLight representative traces 
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/dLight")
for (n in 1:length(csv_files)) {
  
  plot <- ggplot(Filtered_Traces[[n]], aes(x=Time, y=FoldChange)) +
    geom_line() +
    xlim(0,1500) + xlab("Time (sec.)") + scale_x_continuous(breaks = seq(0, 1500, by = 300)) +
    ylim(0.8, 1.2) + 
    ylab("dLight Fluorescence (AU)") + 
    theme_classic()
  
  ggsave(paste(csv_files[n], "_Trace.pdf", sep = ""),
         width = 8, 
         height = 4, 
         plot = plot)
}

#Export the dLight representative traces for 100 sec.
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/dLight")
for (n in 1:length(csv_files)) {
  
  plot <- ggplot(Filtered_Traces[[n]], aes(x=Time, y=FoldChange)) +
    geom_line() +
    xlim(700,800) + 
    xlab("Time (sec.)") + 
    #scale_x_continuous(breaks = 20) +
    ylim(0.8, 1.2) + 
    ylab("dLight Fluorescence (AU)") + 
    theme_classic()
  
  ggsave(paste(csv_files[n], "_100sec_Trace.pdf", sep = ""),
         width = 8, 
         height = 4, 
         plot = plot)
}
#-------------------------------------------------------------------------------
#Normalize the curve based on the fold change of the first bin
max_time <- 1500
bin_width <- 100

#Bin the DeltaF/F values
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]]$Bin <- cut(Filtered_Traces[[n]]$Time, breaks = seq(from = 0, to = max_time, by = bin_width)) #1500/30 = 50 bins
  Filtered_Traces[[n]]$Bin <- as.numeric(Filtered_Traces[[n]]$Bin)
}

#Generate a dataframe to store the FC per bin for each mouse
dLight_FC <- data.frame(matrix(NA, nrow = max_time/bin_width, ncol = length(csv_files)))
rownames(dLight_FC) <- unique(rownames(Filtered_Traces[[1]]$Bin))
colnames(dLight_FC) <- csv_files

#Compute the FC of dLight signal based on the first bin
for (Bin_x in 1:nrow(dLight_FC)) {
  for (sample_y in 1:ncol(dLight_FC)) {
    dLight_FC[Bin_x, sample_y] <- mean(Filtered_Traces[[sample_y]][Filtered_Traces[[sample_y]]$Bin == Bin_x, ]$AIN01) / mean(Filtered_Traces[[sample_y]][Filtered_Traces[[sample_y]]$Bin == 1, ]$AIN01)
  }
}

#Stack the dataframe
dLight_FC <- stack(dLight_FC)

#Add a Time column representing the actual bin time values
dLight_FC$Time <- rep(seq(from = bin_width, to = max_time, by = bin_width), times = length(csv_files))

#Add a group column
dLight_FC <- dLight_FC %>%
  mutate(
    Group = case_when(
      ind == "J0782_dLight_0001.csv" ~ "fDIO-mCherry",
      ind == "J0859_dLight_0001.csv" ~ "fDIO-mCherry",
      ind == "J0777_dLight_0001.csv" ~ "fDIO-mCherry",
      ind == "K0992_dLight_0001.csv" ~ "fDIO-mCherry",
      ind == "K0871_dLight_0001.csv" ~ "fDIO-mCherry",
      
      ind == "J0784_dLight_0001.csv" ~ "CoffFon-DREADD",
      ind == "J0861_dLight_0001.csv" ~ "CoffFon-DREADD",
      ind == "J0779_dLight_0001.csv" ~ "CoffFon-DREADD",
      ind == "K0331_dLight_0001.csv" ~ "CoffFon-DREADD",
      ind == "K0994_dLight_0001.csv" ~ "CoffFon-DREADD",
      
      ind == "K0023_dLight_0001.csv" ~ "ConFon-DREADD",
      ind == "K0029_dLight_0001.csv" ~ "ConFon-DREADD",
      ind == "K0863_dLight_0001.csv" ~ "ConFon-DREADD",
      ind == "K0865_dLight_0001.csv" ~ "ConFon-DREADD",
      ind == "K0990_dLight_0001.csv" ~ "ConFon-DREADD",
      
    )
  )

dLight_FC |>
  tidyplot(x = Time, y = values, color = Group, dodge_width = 0) |>
  add_mean_line() |>
  add_ci95_ribbon() |>
  adjust_x_axis_title("Time (sec.)") |>
  adjust_x_axis(breaks = seq(from = 0, to = 1500, by = 300)) |>
  adjust_y_axis_title("dLight Fold Change") |>
  adjust_y_axis(limits = c(0.82,1.07), breaks = seq(from = 0.8, to = 1.1, by = 0.05))

#-------------------------------------------------------------------------------
################################################################################
##########################  PROCESSING OF ISO TRACES  ##########################
################################################################################
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/Iso")

#Import the traces as a list
csv_files <- list.files(pattern = "\\.csv$")

Traces <- list()
for (n in 1:length(csv_files)) {
  Traces[[n]] <- read.csv(csv_files[n])
}

names(Traces) <- csv_files

#-------------------------------------------------------------------------------
#Export the traces as graph
#Filter the traces to be included only in the first 25 min
Filtered_Traces <- Traces
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]] <- Traces[[n]][Traces[[n]]$Time < 1500 & Traces[[n]]$Time > 5,]
}

#Normalise the curve value on the average value of the first 100 seconds
for (mouseID in names(Filtered_Traces)) {
  Filtered_Traces[[mouseID]]$FoldChange <- Filtered_Traces[[mouseID]]$AIN01 /mean(Filtered_Traces[[mouseID]][Filtered_Traces[[mouseID]]$Time < 100,]$AIN01)
  
}

#Export the DDF dLight representative traces 
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/Iso")
for (n in 1:length(csv_files)) {
  
  plot <- ggplot(Filtered_Traces[[n]], aes(x=Time, y=FoldChange)) +
    geom_line() +
    xlim(0,1500) + xlab("Time (sec.)") + scale_x_continuous(breaks = seq(0, 1500, by = 300)) +
    ylim(0.8, 1.2) + 
    ylab("dLight Fluorescence (AU)") + 
    theme_classic()
  
  ggsave(paste(csv_files[n], "_Trace.pdf", sep = ""),
         width = 8, 
         height = 4, 
         plot = plot)
}

#Export the iso representative traces for 100 sec.
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/Iso")
for (n in 1:length(csv_files)) {
  
  plot <- ggplot(Filtered_Traces[[n]], aes(x=Time, y=FoldChange)) +
    geom_line() +
    xlim(700,800) + 
    xlab("Time (sec.)") + 
    #scale_x_continuous(breaks = 20) +
    ylim(0.8, 1.2) + 
    ylab("dLight Fluorescence (AU)") + 
    theme_classic()
  
  ggsave(paste(csv_files[n], "_100sec_Trace.pdf", sep = ""),
         width = 8, 
         height = 4, 
         plot = plot)
}

#-------------------------------------------------------------------------------
#Normalize the curve based on the fold change of the first bin
max_time <- 1500
bin_width <- 100

#Bin the DeltaF/F values
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]]$Bin <- cut(Filtered_Traces[[n]]$Time, breaks = seq(from = 0, to = max_time, by = bin_width)) #1500/30 = 50 bins
  Filtered_Traces[[n]]$Bin <- as.numeric(Filtered_Traces[[n]]$Bin)
}

#Generate a dataframe to store the FC per bin for each mouse
Iso_FC <- data.frame(matrix(NA, nrow = max_time/bin_width, ncol = length(csv_files)))
rownames(Iso_FC) <- unique(rownames(Filtered_Traces[[1]]$Bin))
colnames(Iso_FC) <- csv_files

#Compute the FC of Iso signal based on the first bin
for (Bin_x in 1:nrow(Iso_FC)) {
  for (sample_y in 1:ncol(Iso_FC)) {
    Iso_FC[Bin_x, sample_y] <- mean(Filtered_Traces[[sample_y]][Filtered_Traces[[sample_y]]$Bin == Bin_x, ]$AIN01) / mean(Filtered_Traces[[sample_y]][Filtered_Traces[[sample_y]]$Bin == 1, ]$AIN01)
  }
}

#Stack the dataframe
Iso_FC <- stack(Iso_FC)

#Add a Time column representing the actual bin time values
Iso_FC$Time <- rep(seq(from = bin_width, to = max_time, by = bin_width), times = length(csv_files))

#Add a group column
Iso_FC <- Iso_FC %>%
  mutate(
    Group = case_when(
      ind == "J0782_Iso_0000.csv" ~ "fDIO-mCherry",
      ind == "J0859_Iso_0000.csv" ~ "fDIO-mCherry",
      ind == "J0777_Iso_0000.csv" ~ "fDIO-mCherry",
      ind == "K0992_Iso_0000.csv" ~ "fDIO-mCherry",
      ind == "K0871_Iso_0000.csv" ~ "fDIO-mCherry",
      
      ind == "J0784_Iso_0000.csv" ~ "CoffFon-DREADD",
      ind == "J0861_Iso_0000.csv" ~ "CoffFon-DREADD",
      ind == "J0779_Iso_0000.csv" ~ "CoffFon-DREADD",
      ind == "K0331_Iso_0000.csv" ~ "CoffFon-DREADD",
      ind == "K0994_Iso_0000.csv" ~ "CoffFon-DREADD",
      
      ind == "K0023_Iso_0000.csv" ~ "ConFon-DREADD",
      ind == "K0029_Iso_0000.csv" ~ "ConFon-DREADD",
      ind == "K0863_Iso_0000.csv" ~ "ConFon-DREADD",
      ind == "K0865_Iso_0000.csv" ~ "ConFon-DREADD",
      ind == "K0990_Iso_0000.csv" ~ "ConFon-DREADD",
      
    )
  )

Iso_FC |>
  tidyplot(x = Time, y = values, color = Group, dodge_width = 0) |>
  add_mean_line() |>
  add_ci95_ribbon() |>
  adjust_x_axis_title("Time (sec.)") |>
  adjust_x_axis(breaks = seq(from = 0, to = 1500, by = 300)) |>
  adjust_y_axis_title("Isobestic Fold Change") |>
  adjust_y_axis(limits = c(0.82,1.07), breaks = seq(from = 0.8, to = 1.1, by = 0.05))



