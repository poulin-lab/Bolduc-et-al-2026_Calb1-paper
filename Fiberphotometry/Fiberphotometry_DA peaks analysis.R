library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(tidyplots)
library(purrr)
library(zoo)

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

names(Traces) <- csv_files

#-------------------------------------------------------------------------------
#Export the traces as graph
#Filter the traces to be included only in the first 25 min
Time_Threshold <- 2100

Filtered_Traces <- Traces
for (n in 1:length(csv_files)) {
  Filtered_Traces[[n]] <- Traces[[n]][Traces[[n]]$Time < Time_Threshold & Traces[[n]]$Time > 5,]
}

#Add AIN01 which is the motion corrected DFF (%)
for (mouseID in names(Filtered_Traces)) {
  Filtered_Traces[[mouseID]]$AIN01 <- Filtered_Traces[[mouseID]]$DFF * 100
}

#-------------------------------------------------------------------------------
#IMPORTANT NOTE: DO NOT DO IF EXTRACTING THE PEAKS BASED ON THE LOESS FIT

#Split the DFF curve into 100sec. bins
Bins_length <- 25

for (mouseID in names(Filtered_Traces)) {
  time_vector <- Filtered_Traces[[mouseID]]$Time
  breaks <- seq(0, Time_Threshold, by = Bins_length)
  Filtered_Traces[[mouseID]]$binID <- as.numeric(cut(time_vector, 
                                                     breaks = breaks, 
                                                     include.lowest = TRUE, 
                                                     right = FALSE))
  time_vector <- NULL
  breaks <- NULL
                          
}

#-------------------------------------------------------------------------------
#IMPORTANT NOTE: DO NOT DO IF EXTRACTING THE PEAKS BASED ON THE LOESS FIT

#Extracting the peaks of DA release
#Note: A peak will be defined as a time interval where all deltaF/F values are > Threshold % + mode distribution of 100 sec. deltaF/F

DFF_Perc_Threshold <- 1 #Threshold to set to consider it as peak 
Peaks <- list()

for (mouseID in names(Filtered_Traces)) {
  Peaks[[mouseID]] <- list()
  
  # Get unique bins, excluding NA if any
  bins <- unique(Filtered_Traces[[mouseID]]$binID)
  bins <- bins[!is.na(bins)]
  
  for (Bin_n in bins) {
    
    MouseDFF <- Filtered_Traces[[mouseID]] %>%
      dplyr::filter(binID == Bin_n) # Fixed typo: binID instead of BinID
    
    # Define the threshold using kernel density mode
    density_obj <- density(MouseDFF$AIN01)
    Mode_DFF <- density_obj$x[which.max(density_obj$y)]
    Threshold <- Mode_DFF + DFF_Perc_Threshold
    
    Peaks[[mouseID]][["Bins Information"]][[Bin_n]] <- data.frame(Mode_DFF, Threshold)
    
    # Extract the peaks
    idx <- which(MouseDFF$AIN01 > Threshold)
    
    # CHECK: Only proceed if peaks are found in this bin
    if (length(idx) > 0) {
      grp <- cumsum(c(1, diff(idx) != 1)) 
      list_of_peaks <- split(idx, grp)
      Peaks[[mouseID]][["Peaks"]][[Bin_n]] <- lapply(list_of_peaks, function(i) MouseDFF[i, ])
      
      Start <- numeric(length(list_of_peaks))
      End <- numeric(length(list_of_peaks))
      
      for (n in 1:length(list_of_peaks)) {
        Peak_n <- Peaks[[mouseID]][["Peaks"]][[Bin_n]][[n]]
        Start[n] <- Peak_n$Time[1]
        End[n] <- Peak_n$Time[nrow(Peak_n)]
      }
      
      Peaks[[mouseID]][["PerBin_Start/End"]][[Bin_n]] <- data.frame(Start = Start, End = End)
      
    } else {
      # Handle cases with no peaks to avoid the data.frame error
      Peaks[[mouseID]][["Peaks"]][[Bin_n]] <- list()
      Peaks[[mouseID]][["PerBin_Start/End"]][[Bin_n]] <- data.frame(Start = numeric(0), End = numeric(0))
    }
  }
}

#Add a column for the mode value in every peak
for (mouseID in names(Peaks)) {
  for (Bin_n in 1:length(Peaks[[mouseID]][["Peaks"]])) {
    if (length(Peaks[[mouseID]][["Peaks"]][[Bin_n]]) > 0) {
      for (Peak_n in 1:length(Peaks[[mouseID]][["Peaks"]][[Bin_n]])) {
        Peaks[[mouseID]][["Peaks"]][[Bin_n]][[Peak_n]]$ModeDFF <-  rep(Peaks[[mouseID]][["Bins Information"]][[Bin_n]]$Mode_DFF,
                                                                       times = nrow(Peaks[[mouseID]][["Peaks"]][[Bin_n]][[Peak_n]]))
                                                          
      }
    }
  }
}

#Stack all traces per animal
for (mouseID in names(Peaks)) {
  Peaks[[mouseID]][["All_Peaks_Combined"]] <- purrr::list_flatten(Peaks[[mouseID]][["Peaks"]])
  Peaks[[mouseID]][["Start/End"]] <- dplyr::bind_rows(Peaks[[mouseID]][["PerBin_Start/End"]])
}

#Add the mode per binID to the Filtered_Traces
for (mouseID in names(Filtered_Traces)) {
  
  # 1. Extract and stack the Bin Information for this mouse
  bin_info_df <- dplyr::bind_rows(Peaks[[mouseID]][["Bins Information"]], .id = "binID")
  
  # 2. Ensure binID is numeric so it matches your Filtered_Traces$binID
  bin_info_df$binID <- as.numeric(bin_info_df$binID)
  
  # 3. Join the info back to the main trace data
  Filtered_Traces[[mouseID]] <- Filtered_Traces[[mouseID]] %>%
    dplyr::left_join(bin_info_df, by = "binID")
  
}

#-------------------------------------------------------------------------------
#Extracting the peaks of DA release based on the LOESS
#Note: alpha = 0.1 showed the best fitting based on a preliminary screening

#Perform a Local Polynomial Regression Fitting (LOESS) fit of the DFF trace
for (mouseID in names(Filtered_Traces)) {
  
  print(paste("Performing", mouseID, "loess alpha = 0.1"))
  
  loess_fit <- loess(formula = AIN01 ~ Time, 
                     data = Filtered_Traces[[mouseID]], 
                     model = FALSE,
                     span = 0.1, #the alpha parameter which controls the degree of smoothing
                     degree = 2,
                     parametric = FALSE, 
                     drop.square = FALSE, 
                     normalize = TRUE,
                     family = "symmetric",
                     method = "loess")
  
  Filtered_Traces[[mouseID]]$loess0.1 <- predict(loess_fit)
}

#Note: A peak will be defined as a time interval where all deltaF/F values are > Threshold % + mode distribution of 100 sec. deltaF/F
DFF_Perc_Threshold <- 1 #Threshold to set to consider it as peak 

#Add a threshold column
for (mouseID in names(Filtered_Traces)) {
  Filtered_Traces[[mouseID]]$Threshold <- Filtered_Traces[[mouseID]]$loess0.1 + DFF_Perc_Threshold
} 

#Extract the peaks
Peaks <- list()

for (mouseID in names(Filtered_Traces)) {
  
  # 1. Initialize lists for this mouse
  Peaks[[mouseID]] <- list()
  
  # 2. Identify all time points above the dynamic LOESS threshold
  idx <- which(Filtered_Traces[[mouseID]]$AIN01 > Filtered_Traces[[mouseID]]$Threshold)
  
  if (length(idx) > 0) {
    
    # 3. Group consecutive indices into distinct peaks
    grp <- cumsum(c(1, diff(idx) != 1)) 
    
    # 4. Split the indices into a list of peaks
    peak_list <- split(idx, grp)
    
    # 5. Extract the full data for each peak
    Peaks[[mouseID]][["Peaks_Data"]] <- lapply(peak_list, function(i) Filtered_Traces[[mouseID]][i, ])
    
    # 6. Extract Start/End times and Summary
    Peaks[[mouseID]][["Summary"]] <- purrr::map_df(Peaks[[mouseID]][["Peaks_Data"]], function(df) {
      data.frame(
        Start = min(df$Time),
        End = max(df$Time),
        Duration = max(df$Time) - min(df$Time),
        Max_Amplitude = max(df$AIN01),
        Peak_Height_Above_Threshold = max(df$AIN01 - df$Threshold),
        Area_Under_Curve = sum(df$AIN01 - df$Threshold)
      )
    })
    
    message(paste(mouseID, ": Found", nrow(Peaks[[mouseID]][["Summary"]]), "peaks."))
    
  } else {
    Peaks[[mouseID]][["Peaks_Data"]] <- list()
    Peaks[[mouseID]][["Summary"]] <- data.frame()
    message(paste(mouseID, ": No peaks found."))
  }
}

#-------------------------------------------------------------------------------
#Export the traces
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peaks analysis")

for (mouseID in names(Filtered_Traces)) {
  
  plot <- ggplot(Filtered_Traces[[mouseID]], aes(x=Time)) +
    geom_line(aes(y = AIN01, color = "dLight DFF (%)"), linewidth = 0.5) +
    geom_line(aes(y = loess0.1, color = "loess_alpha=0.1"), linewidth = 1) +
    geom_line(aes(y = Threshold, color = "Threshold"), linetype = "dotted", linewidth = 1) +
    
    
    # 4. Customizing the colors and legend
    scale_color_manual(name = "Traces", 
                       values = c("dLight DFF (%)" = "black", 
                                  "loess_alpha=0.1" = "blue",
                                  "Threshold" = "red"
)) +
    
    scale_x_continuous(limits = c(0, 1500), 
                       breaks = seq(0, 1500, by = 10)) +
    
    xlab("Time (sec.)") + 
    ylim(-15, 15) + 
    ylab("dLight deltaF/F (%)") + 
    theme_classic() +
    geom_rect(data = Peaks[[mouseID]][["Summary"]], mapping = aes(xmin = Start, xmax = End, ymin = -15, ymax = 15), alpha = 0.4, inherit.aes = FALSE) 
  
  ggsave(paste(mouseID, "_PeaksLabelled.pdf", sep = ""),
         width = 150, 
         height = 4, 
         units = "in", 
         limitsize = FALSE,
         plot = plot)
}

for (mouseID in names(Filtered_Traces)) {
  
  plot <- ggplot(Filtered_Traces[[mouseID]], aes(x=Time)) +
    geom_line(aes(y = AIN01, color = "dLight DFF (%)"), linewidth = 0.5) +
    geom_line(aes(y = loess0.1, color = "loess_alpha=0.1"), linewidth = 1) +
    geom_line(aes(y = Threshold, color = "Threshold"), linetype = "dotted", linewidth = 1) +
    
    
    # 4. Customizing the colors and legend
    scale_color_manual(name = "Traces", 
                       values = c("dLight DFF (%)" = "black", 
                                  "loess_alpha=0.1" = "blue",
                                  "Threshold" = "red"
                       )) +
    
    scale_x_continuous(limits = c(10, 110), 
                       breaks = seq(0, 1500, by = 10)) +
    
    xlab("Time (sec.)") + 
    ylim(-15, 15) + 
    ylab("dLight deltaF/F (%)") + 
    theme_classic() +
    geom_rect(data = Peaks[[mouseID]][["Summary"]], mapping = aes(xmin = Start, xmax = End, ymin = -15, ymax = 15), alpha = 0.4, inherit.aes = FALSE) 
  
  ggsave(paste(mouseID, "_10-110sec_PeaksLabelled.pdf", sep = ""),
         width = 25, 
         height = 4, 
         units = "in", 
         limitsize = FALSE,
         plot = plot)
}

for (mouseID in names(Filtered_Traces)) {
  
  plot <- ggplot(Filtered_Traces[[mouseID]], aes(x=Time)) +
    geom_line(aes(y = AIN01, color = "dLight DFF (%)"), linewidth = 0.5) +
    geom_line(aes(y = loess0.1, color = "loess_alpha=0.1"), linewidth = 1) +
    geom_line(aes(y = Threshold, color = "Threshold"), linetype = "dotted", linewidth = 1) +
    
    
    # 4. Customizing the colors and legend
    scale_color_manual(name = "Traces", 
                       values = c("dLight DFF (%)" = "black", 
                                  "loess_alpha=0.1" = "blue",
                                  "Threshold" = "red"
                       )) +
    
    scale_x_continuous(limits = c(1000, 1100), 
                       breaks = seq(0, 1500, by = 10)) +
    
    xlab("Time (sec.)") + 
    ylim(-15, 15) + 
    ylab("dLight deltaF/F (%)") + 
    theme_classic() +
    geom_rect(data = Peaks[[mouseID]][["Summary"]], mapping = aes(xmin = Start, xmax = End, ymin = -15, ymax = 15), alpha = 0.4, inherit.aes = FALSE) 
  
  ggsave(paste(mouseID, "_1000-1100sec_PeaksLabelled.pdf", sep = ""),
         width = 25, 
         height = 4, 
         units = "in", 
         limitsize = FALSE,
         plot = plot)
}

#-------------------------------------------------------------------------------
#Plot the peaks amplitude frequency pre and post-dcz injection
Plot_PeakFreq <- function(Start_Time, End_Time) {
  
  Time_Range <- End_Time - Start_Time
  
  plot_data_list <- list()
  
  for (mouseID in names(Peaks)) {
    #Extract the peaks found within the time frame
    df <- Peaks[[mouseID]][["Summary"]] %>%
      dplyr::filter(Start >= Start_Time & Start <= End_Time)
    
    #Extract the frequency
    breaks <- seq(0, 22, by = 0.25)
    Peak_Amplitude <- df$Peak_Height_Above_Threshold + DFF_Perc_Threshold
    intervals <- cut(Peak_Amplitude, breaks = breaks, right = FALSE)
    counts <- as.data.frame(table(intervals))
    counts$Freq <- counts$Freq/Time_Range 
    counts$Amplitude <- breaks[1:(length(breaks)-1)]
    counts$mouseID <- rep(unlist(strsplit(mouseID, split = "_DFF_0000.csv")), times = nrow(counts))
    plot_data_list[[mouseID]] <- counts
    
  }
  
  Toplot <- dplyr::bind_rows(plot_data_list)
  
  #Add a group column
  Toplot <- Toplot %>%
    mutate(
      Group = case_when(
        mouseID == "J0782" ~ "fDIO-mCherry",
        mouseID == "J0859" ~ "fDIO-mCherry",
        mouseID == "J0777" ~ "fDIO-mCherry",
        mouseID == "K0992" ~ "fDIO-mCherry",
        mouseID == "K0871" ~ "fDIO-mCherry",
        mouseID == "J0871" ~ "fDIO-mCherry",
        
        mouseID == "J0784" ~ "CoffFon-DREADD",
        mouseID == "J0861" ~ "CoffFon-DREADD",
        mouseID == "K0861" ~ "CoffFon-DREADD",
        mouseID == "J0779" ~ "CoffFon-DREADD",
        mouseID == "K0331" ~ "CoffFon-DREADD",
        mouseID == "K0994" ~ "CoffFon-DREADD",
        
        mouseID == "K0023" ~ "ConFon-DREADD",
        mouseID == "K0029" ~ "ConFon-DREADD",
        mouseID == "K0863" ~ "ConFon-DREADD",
        mouseID == "K0865" ~ "ConFon-DREADD",
        mouseID == "K0990" ~ "ConFon-DREADD",
        
      )
    )
  
  
  #Plot the graph of dLight DFF SD
  Toplot |>
    tidyplot(x = Amplitude, y = Freq, color = Group, dodge_width = 0) |>
    add_mean_line() |>
    add_ci95_ribbon() |>
    adjust_x_axis_title("Peak Amplitude (%)") |>
    adjust_x_axis(limits = c(DFF_Perc_Threshold, 10), breaks = seq(from = 0, to = 15, by = 2.5)) |>
    adjust_y_axis(limits = c(0, 0.17), breaks = seq(from = 0, to = 0.15, by = 0.05)) |>
    adjust_y_axis_title("Frequency (Hz)")
  
}

Plot_PeakFreq(Start_Time = 0,
              End_Time = 300)

Plot_PeakFreq(Start_Time = 500,
              End_Time = 1500)

write_rds(Filtered_Traces, "LOESS_Traces.rds")


