library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(tidyplots)
library(MASS)
library(signal)
library(viridis)

#----------------------------------------------------------------
dLight_Directory <- "C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/dLight"

Iso_Directory <- "C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/Iso"

#---------------------------------------------------------------
#Read the csv
#Import the traces as a list
Load_Traces <- function(Directory) {
  
  setwd(Directory)
  
  csv_files <- list.files(pattern = "\\.csv$")
  Traces <- list()
  for (n in 1:length(csv_files)) {
    Traces[[n]] <- read.csv(csv_files[n])
  }
  
  csv_files <- sapply(strsplit(csv_files, split =  "_"), "[[", 1)
  names(Traces) <- csv_files
  return(Traces)

}

dLight_Traces <- Load_Traces(Directory = dLight_Directory)
Iso_Traces <- Load_Traces(Directory = Iso_Directory)

#Combine into one list
Traces <- dLight_Traces
for (mouseID in names(dLight_Traces)) {
  colnames(Traces[[mouseID]]) <- c("f470", "Time")
  Traces[[mouseID]]$f405 <- Iso_Traces[[mouseID]]$AIN01
  Traces[[mouseID]] <- dplyr::filter(Traces[[mouseID]], Time < 2100) #Filter the times that is lower than 2100
}

dLight_Traces <- NULL
Iso_Traces <- NULL

#---------------------------------------------------------------
#Export the traces

for (mouseID in names(Traces)) {
  setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/dLight")
  
  Toplot <- Traces[[mouseID]]
  
plot <- ggplot(data = Toplot, mapping = aes(x = Time, y = f470)) +
    geom_line() +
    xlab("Time (s)") +
    ylab("dLight F") + 
    theme_classic()
  
  ggsave(
    paste(mouseID, "dLightTrace.pdf", sep = ""),
    plot = plot,
    width = 20, 
    height = 5, 
    units = "cm")
    
  
  
    setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/Iso")
    
    Toplot <- Traces[[mouseID]]
    
    plot <- ggplot(data = Toplot, mapping = aes(x = Time, y = f405)) +
      geom_line() +
      xlab("Time (s)") +
      ylab("dLight F") + 
      theme_classic()
    
    ggsave(
      paste(mouseID, "IsoTrace.pdf", sep = ""),
      plot = plot,
      width = 20, 
      height = 5, 
      units = "cm")
}

#---------------------------------------------------------------
#Remove the events that represent artifacts

#Disconnections occured at these time for K0331
Traces[["K0331"]] <- Traces[["K0331"]] |> 
  dplyr::filter(!between(Time, 1640, 1680))

Traces[["K0331"]] <- Traces[["K0331"]] |> 
  dplyr::filter(!between(Time, 1820, 1860))

ggplot(data = Toplot, mapping = aes(x = Time, y = f405)) +
  geom_line() +
  xlab("Time (s)") +
  ylab("dLight F") + 
  theme_classic()

#---------------------------------------------------------------
###############################################################
################## Motion correction ##########################
###############################################################
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Traces/DFF/LowBF3Hz_IRLS300sec")

for (mouseID in names(Traces)) {
  Data <- Traces[[mouseID]]
  Time <- Data$Time
  
  fs <- 60 # fs: Frequency of sampling (Hz) i.e nFrames per second
  f405 <- Data$f405 #vector of isosbestic signal
  f470 <- Data$f470 #vector of dLight signal
  
  # 1. Low-pass butterworth filter
  bf <- butter(2, 3 / (fs/2), type = "low") # Low-pass butterworth filter of 2nd degree, 3Hz
  f470_filt <- filtfilt(bf, f470)
  f405_filt <- filtfilt(bf, f405)
  
  # 2. IRLS regresion to baseline
  df_baseline <- data.frame(f470_filt = f470_filt[Time > 5 & Time < 300], 
                            f405_filt = f405_filt[Time > 5 & Time < 300])
  
  fit_baseline <- rlm(f470_filt ~ f405_filt, data = df_baseline, psi = psi.huber, k = 1.4)
  f405_fitted <- predict(fit_baseline, newdata = data.frame(f405_filt = f405_filt))
  
  # 3. deltaF/F calculation
  DFF <- (f470_filt - f405_fitted) / f405_fitted
  
  Traces[[mouseID]]$f470_filt <- f470_filt
  Traces[[mouseID]]$f405_fitted <- f405_fitted
  Traces[[mouseID]]$DFF <- DFF
  
  # 4. Visualisation
  Toplot <- data.frame(Time = Time,
                       f405 = f405,
                       f470 = f470,
                       f405_filt = f405_filt,
                       f470_filt = f470_filt,
                       f405_fitted = f405_fitted,
                       DFF = DFF)
  
  # Frame the data frame to the correct orientation for ggplot
  Toplot <- Toplot %>%
    pivot_longer(
      cols = -Time,           # Keep the Time column to pivot the rest
      names_to = "Identity",   # Name of the identity column
      values_to = "Value"      # Name of the numeric column
    )
  
  Toplot <- Toplot %>%
    dplyr::filter(Identity %in% c(
      #"f405",
      #"f470",
      #"f405_filt",
      "f470_filt",
      "f405_fitted",
      "DFF"
    ))
  
  Toplot$Identity <- factor(Toplot$Identity, 
                            levels = c(
                              #"f405",
                              #"f470",
                              #"f405_filt",
                              "DFF",
                              "f470_filt",
                              "f405_fitted"
                              
                            ))
  
  plot <- ggplot(Toplot, mapping = aes(x = Time, y = Value, color = Identity)) +
    xlim(0, 2100) +
    ylim(-0.2, 0.2) +
    geom_line(linewidth = 0.1) +
    theme_classic()  +
    scale_color_viridis_d(option = "turbo")
  
  ggsave(
    paste(mouseID, "DFFTrace.pdf", sep = ""),
    plot = plot,
    width = 20, 
    height = 5, 
    units = "cm")
  
}


#---------------------------------------------------------------
#Export the traces
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/DFF/LowBF3Hz_IRLS300sec")
for (mouseID in names(Traces)) {
  write.csv(Traces[[mouseID]], paste(mouseID, "_DFF_0000.csv", sep = ""))
}

