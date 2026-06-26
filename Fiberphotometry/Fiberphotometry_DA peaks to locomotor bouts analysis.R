library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(tidyplots)
library(data.table)
library(DescTools)
library(magick)
library(imager)
library(lattice)
library(viridis)
#library(paletteer)
library(gplots)
library(sommer)
library(patchwork)
library(gganimate)
library(lubridate)
library(ggpointdensity)
library(stringr)
library(scales)
library(purrr)
library(gglinedensity)
library(ggExtra)
library(WRS2)
library(rstatix)
#-------------------------------------------------------------------------------
################################################################################
##########################  LOADING BODY MOTION TRACES  ########################
################################################################################

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/BodyTraces_csv")

#Import the body traces as a list
csv_files <- list.files(pattern = "\\.csv$")

Body_motion <- list()
for (n in csv_files) {
  Body_motion[[n]] <- read.csv(n)
  Body_motion[[n]]$Time <- Body_motion[[n]]$frame_idx * 0.0333333333333333 #Note add a column to provide the timescale to fibrephotometry, resolution is 0,0333 sec. per frame
}

names(Body_motion) <- unlist(strsplit(csv_files, split = ".analysis.csv"))

#-------------------------------------------------------------------------------
#Load the image
#image_filenames <- list.files(pattern="*.png")
#img <- readPNG(image_filenames[1])
#h <- nrow(img)
#w <- ncol(img)
#windows(width = 10, height = 10 * (h/w)) 
#par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
#plot(1, type = "n", xlim = c(0, w), ylim = c(h, 0), asp = 1, axes = FALSE)
#rasterImage(img, 0, h, w, 0)

#Place two dots in corners
#pts <- locator(2)

#Show the points
#points(pts$x, pts$y, col = "red", pch = 19, cex = 1.5)

#Extract the distance in pixels
#x1 <- pts$x[1]
#y1 <- pts$y[1]
#x2 <- pts$x[2]
#y2 <- pts$y[2]

#Pixel_Distance <- sqrt((x2-x1)^2+(y2-y1)^2)
#Resolution <- 39.5/Pixel_Distance
#Resolution

#Note: The value I got for one pixel is equivalent to 0.04483225, or 0.448mm
Resolution <- 0.448

#-------------------------------------------------------------------------------
#Quality control
Torso <- list()
for (mouseID in names(Body_motion)) {
  Torso[[mouseID]] <- Body_motion[[mouseID]]
  
  #Remove the rows has unlablled torso points
  Torso[[mouseID]] <- Torso[[mouseID]][!is.na(Torso[[mouseID]]$torso.x),]
  Torso[[mouseID]] <- Torso[[mouseID]][!is.na(Torso[[mouseID]]$torso.y),]
  
  #Remove the rows that has no torso score
  Torso[[mouseID]] <- Torso[[mouseID]][!is.na(Torso[[mouseID]]$torso.score),]
  
  #Remove the torso points tnat has a confidence score < 0.8
  Torso[[mouseID]] <- Torso[[mouseID]][Torso[[mouseID]]$torso.score>0.8,]
  
  #Convert the x y coordinates into mm scale
  Torso[[mouseID]][,c("head.x", "head.y", "torso.x", "torso.y", "tail_base.x", "tail_base.y")] <- Torso[[mouseID]][,c("head.x", "head.y", "torso.x", "torso.y", "tail_base.x", "tail_base.y")] * Resolution
  
  #Remove the track column
  Torso[[mouseID]]$track <- NULL
}

#------------------------------------------------------------------------
#Downsample the data to the time resolution
time_resolution <- 0.1 #Time resolution in second

for (mouseID in names(Torso)) {
  print(mouseID)
  Ref_frames <- seq(from = 0, to = 2100, by = time_resolution)
  Closest_frames <- Closest(Torso[[mouseID]]$Time, Ref_frames, which = FALSE, na.rm = FALSE)
  
  New_frames <- c()
  for (x in 1:length(Closest_frames)) {
    New_frames[x] <- Closest_frames[[x]]
  }
  
  New_frames <- unique(New_frames) #Some frames have been assigned in double, so need to take only the unique ones
  
  Frame_indices <- c()
  for (x in 1:length(New_frames)) {
    Frame_indices[x] <- which(Torso[[mouseID]]$Time == New_frames[x])
  }
  
  #Subsample the frames
  Torso[[mouseID]] <- Torso[[mouseID]][Frame_indices,]
  
}

#-------------------------------------------------------------------------------
#Compute the values to velocity
velocity <- function(X1, X2, Y1, Y2, T1, T2) { 
  speed <- (sqrt((X2 - X1)^2 + (Y2-Y1)^2))/(T2-T1)
  return (speed)
}

for (mouseID in names(Torso)) {
  
  print(mouseID)
  Torso.velocity <- c()
  
  for (x in 2:nrow(Torso[[mouseID]])) {
    Torso.velocity[x] <- velocity(X1 = Torso[[mouseID]]$torso.x[x-1],
                                  X2 = Torso[[mouseID]]$torso.x[x],
                                  Y1 = Torso[[mouseID]]$torso.y[x-1],
                                  Y2 = Torso[[mouseID]]$torso.y[x],
                                  T1 = Torso[[mouseID]]$Time[x-1],
                                  T2 = Torso[[mouseID]]$Time[x])
  }
  
  Torso[[mouseID]]$Torso.velocity <- Torso.velocity
  
}

#------------------------------------------------------------------------
#Rearrange the torso list keeping the sleap_data as a separate element
Torso_New <- list()

for (mouseID in names(Torso)) {
  Torso_New[[mouseID]][["Sleap_data"]] <- Torso[[mouseID]]
}

Torso <- Torso_New
Torso_New <- NULL

#------------------------------------------------------------------------
# Extract the bouts
#Note: Bouts will be defined as speed is defined as an occurence where the speed is suddenly increased >25mm/s for at least 600ms 

#Set the threshold for speed and time
Speed_Threshold <- 25 #in mm/s
Time_Threshold <- 0.6 #in sec.

#Extract the bouts

for (mouseID in names(Torso)) {
  
  #Identify the index having a velocity value higher than the time threshold
  idx <- which(Torso[[mouseID]][["Sleap_data"]]$Torso.velocity > Speed_Threshold)
  grp <- cumsum(c(1, diff(idx) != 1)) #Identify index that are seperated by more than one value
  Potential_bouts <- split(idx, grp)
  
  Bouts_Start_Index <- c()
  Bouts_End_Index <- c()
  
  for (n in 1:length(Potential_bouts)) {
    
    Start_Index <- Potential_bouts[[n]][1]
    End_Index <- Potential_bouts[[n]][length(Potential_bouts[[n]])]
    
    #Filter for the bouts that are longer than the time threshold
    if (Torso[[mouseID]][["Sleap_data"]]$Time[End_Index] - Torso[[mouseID]][["Sleap_data"]]$Time[Start_Index] > Time_Threshold) {
      Bouts_Start_Index[n] <- Start_Index - 1 #Adding one index prior occurence which has a speed > threshold
      Bouts_End_Index[n] <- End_Index + 1 #Adding one index post bout occurence which has a speed > threshold
    }
  }
  
  #Remove the NA generated
  Bouts_Start_Index <- Bouts_Start_Index[!is.na(Bouts_Start_Index)]
  Bouts_End_Index <- Bouts_End_Index[!is.na(Bouts_End_Index)]
  
  #Store the bouts index
  Torso[[mouseID]][["Bouts_Index"]] <- data.frame(Start_Index = Bouts_Start_Index, End_Index = Bouts_End_Index)
  
  #Generate a list to store the bouts
  Torso[[mouseID]][["Bouts"]] <- list()
  
  for (Bout_n in 1:nrow(Torso[[mouseID]][["Bouts_Index"]])) {
    Torso[[mouseID]][["Bouts"]][[Bout_n]] <- Torso[[mouseID]][["Sleap_data"]][Torso[[mouseID]][["Bouts_Index"]]$Start_Index[Bout_n]:Torso[[mouseID]][["Bouts_Index"]]$End_Index[Bout_n],]
  }
}

#------------------------------------------------------------------------
#Save the RDS file
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Body_motion")
saveRDS(Torso, file = "Body_motion.rds")

#------------------------------------------------------------------------
#Integrate the body motion and DA traces together
Body_motion <- readRDS("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Body_motion/Body_motion.rds")

#Filter the speed of K0331 at the moment where the cable got unpluged
Body_motion[["K0331"]][["Sleap_data"]] <- Body_motion[["K0331"]][["Sleap_data"]] |> 
  dplyr::filter(!between(Time, 1640, 1680))

Body_motion[["K0331"]][["Sleap_data"]] <- Body_motion[["K0331"]][["Sleap_data"]] |> 
  dplyr::filter(!between(Time, 1820, 1860))

intervals <- list(c(1640, 1680), c(1820, 1860))

Body_motion[["K0331"]][["Bouts"]] <- Filter(function(bout) {
  bout_time <- bout$Time
  
  is_bad <- any(sapply(intervals, function(inter) {
    any(bout_time >= inter[1] & bout_time <= inter[2])
  }))
  
  return(!is_bad)
}, Body_motion[["K0331"]][["Bouts"]])


#Import the traces as a list

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/dLightTraces/DFF/LowBF3Hz_IRLS300sec")
#csv_files <- list.files(pattern = "\\.csv$")
#Traces <- list()
#for (n in 1:length(csv_files)) {
#  Traces[[n]] <- read.csv(csv_files[n])
#}
#names(Traces) <- unlist(strsplit(csv_files, split = "_DFF_0000.csv"))

Traces <- readRDS("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peaks analysis/LOESS_Traces.rds")
names(Traces) <- unlist(strsplit(names(Traces), split = "_DFF_0000.csv"))

#Add the DA trace per mouse
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["dLight_DFF"]] <- Traces[[mouseID]]
}

#There is a mistake of labelling to one of these ID, I will correct it later
#Body_motion[["K0861"]][["dLight_DFF"]] <- Traces[["J0861"]]
Body_motion[["K0871"]][["dLight_DFF"]] <- Traces[["J0871"]]

#Add a AIN01 column that correspond to DFF
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["dLight_DFF"]]$AIN01 <- Body_motion[[mouseID]][["dLight_DFF"]]$DFF * 100
}

#Generate a column where the average trace generated by loess is substracted from the deltaF/F
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["dLight_DFF"]]$DFF_LOESSsub <- Body_motion[[mouseID]][["dLight_DFF"]]$AIN01 - Body_motion[[mouseID]][["dLight_DFF"]]$loess0.1
}

#------------------------------------------------------------------------
#Generate an average trace for dLight prior and post-bout occurence and bout peak-speed

Plot_PeakSpeed <- function(curve_to_timeloc, n_Seconds, min_time, max_time, title, ylim) {
  curve_col_index <- which(colnames(Body_motion[[mouseID]][["dLight_DFF"]]) == curve_to_timeloc)
  Bouts <- list()
  
  for (mouseID in names(Body_motion)) {
    
    Bouts[[mouseID]] <- list()
    Bouts[[mouseID]][["Bouts"]] <- list()
    
    for (Bout_n in seq_along(Body_motion[[mouseID]][["Bouts"]])) {
      start_time <- Body_motion[[mouseID]][["Bouts"]][[Bout_n]]$Time[1]
      
      if (start_time > min_time && start_time < max_time) {
        Bouts[[mouseID]][["Bouts"]][[length(Bouts[[mouseID]][["Bouts"]]) + 1]] <- Body_motion[[mouseID]][["Bouts"]][[Bout_n]]
      }
    }
  }
  
  
  #Time loc these bouts few seconds before and after to their occurence
  #n_Seconds <- 2
  
  for (mouseID in names(Bouts)) {
    
    print(paste("Time loc:", mouseID))
    
    
    Bouts[[mouseID]][["Time_Range"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      Bout <- Bouts[[mouseID]][["Bouts"]][[Bout_n]]
      Start_Time <- Bout$Time[1] 
      Pre_Time <- Start_Time - n_Seconds 
      Post_Time <- Start_Time + n_Seconds
      Peak_Time <- Bout[which.max(Bout$Time),]$Time
      Peak_Pre_Time <- Peak_Time - n_Seconds
      Peak_Post_Time <- Peak_Time + n_Seconds
      
      Bouts[[mouseID]][["Time_Range"]][[Bout_n]] <- data.frame(Start_Time = Start_Time, 
                                                               Pre_Time = Pre_Time, 
                                                               Post_Time = Post_Time,
                                                               Peak_Time = Peak_Time,
                                                               Peak_Pre_Time = Peak_Pre_Time,
                                                               Peak_Post_Time = Peak_Post_Time) 
    }
  }
  
  
  #Time loc the dLight traces
  for (mouseID in names(Bouts)) {
    
    print(paste("Time loc step2:", mouseID))
    
    Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]] <- list()
    Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      Start_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Start_Time
      Pre_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Pre_Time
      Post_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Post_Time
      Peak_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Time
      Peak_Pre_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Pre_Time
      Peak_Post_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Post_Time
      
      dLight_DFF <- Body_motion[[mouseID]][["dLight_DFF"]]
      dLight_DFF <- dLight_DFF[dLight_DFF$Time > Pre_Time & dLight_DFF$Time < Post_Time,]
      dLight_DFF$Time_Loc <-  dLight_DFF$Time - Start_Time
      Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]][[Bout_n]] <- dLight_DFF
      
      dLight_DFF <- Body_motion[[mouseID]][["dLight_DFF"]]
      dLight_DFF <- dLight_DFF[dLight_DFF$Time > Peak_Pre_Time & dLight_DFF$Time < Peak_Post_Time,]
      dLight_DFF$Time_Loc <-  dLight_DFF$Time - Peak_Time
      Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]][[Bout_n]] <- dLight_DFF
      
    }
  }
  
  #Since the T times for Time loc are not exactly the same, interpolate the deltaF/F values at consensus times
  for (mouseID in names(Bouts)) {
    
    print(paste("Interpolation:", mouseID))
    

    Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]] <- list()
    Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {

      Time <- seq(from = -n_Seconds, to = n_Seconds, by = 0.02)
      data <- Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]][[Bout_n]]
      Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]] <- as.data.frame(approx(x = data$Time_Loc, y = data[,curve_col_index], xout = Time))
      colnames(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]]) <- c("Time", "dLight_DFF")
      
      Time <- seq(from = -n_Seconds, to = n_Seconds, by = 0.02)
      data <- Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]
      Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]] <- as.data.frame(approx(x = data$Time_Loc, y = data[,curve_col_index], xout = Time))
      colnames(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]) <- c("Time", "dLight_DFF")
    }
  }
  
  #Combine all the interpolated traces
  All_traces <- list()
  counter <- 1
  for (mouseID in names(Bouts)) {
    for (Bout_n in seq_along(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]])) {
      All_traces[[counter]] <- Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]]
      All_traces[[counter]]$TraceID <- as.character(rep(counter, times = nrow(All_traces[[counter]])))
      All_traces[[counter]]$mouseID <- rep(mouseID, times = nrow(All_traces[[counter]]))
      counter <- counter + 1
    }
  }
  
  from_BoutStart <- bind_rows(All_traces)
  from_BoutStart <- from_BoutStart %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  All_traces <- list()
  counter <- 1
  for (mouseID in names(Bouts)) {
    for (Bout_n in seq_along(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]])) {
      All_traces[[counter]] <- Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]
      All_traces[[counter]]$TraceID <- as.character(rep(counter, times = nrow(All_traces[[counter]])))
      All_traces[[counter]]$mouseID <- rep(mouseID, times = nrow(All_traces[[counter]]))
      counter <- counter + 1
    }
  }
  
  from_PeakSpeed <- bind_rows(All_traces)
  from_PeakSpeed <- from_PeakSpeed %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )

  plot <- from_PeakSpeed |>
    tidyplot(x = Time, y = dLight_DFF, color = Group) |>
    add_mean_line() |>
    add_ci95_ribbon() |>
    adjust_title(title) |>
    adjust_x_axis_title("Time from peak speed (sec.)") |>
    adjust_x_axis(limits = c(-2, 1.5), breaks = seq(from = -10, to = 10, by = 0.5)) |>
    adjust_y_axis(limits = ylim, breaks = seq(from = -10, to = 10, by = 1)) |>
    adjust_y_axis_title("dLight dF/F (%)") #|>
  
  ggsave(paste(title, "_PeakSpeed.pdf", sep = ""),
         width = 4, 
         height = 4, 
         plot = plot)

  plot
  
  #Collect the number of bouts per mouse
  nBouts <- c()
  for (mouseID_n in 1:length(Bouts)) {
    nBouts[mouseID_n] <- length(Bouts[[mouseID_n]][["Bouts"]])
  }
  
  nBouts_df <- data.frame(mouseID = names(Bouts), nBouts = nBouts)
  nBouts_df <- nBouts_df %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  print(nBouts_df)
  write.csv(nBouts_df, paste(title, "nBouts.csv", sep = "_")) #save the data to CSV
  
  nBouts_df_pergroup <- aggregate(nBouts ~ Group, data = nBouts_df, FUN = sum, na.rm = TRUE)
  print(nBouts_df_pergroup)
  write.csv(nBouts_df_pergroup, paste(title, "nBoutsperGroup.csv", sep = "_")) #save the data to CSV
  

}

Plot_BoutStart <- function(curve_to_timeloc, n_Seconds, ylim, min_time, max_time, title ) {
  curve_col_index <- which(colnames(Body_motion[[mouseID]][["dLight_DFF"]]) == curve_to_timeloc)
  Bouts <- list()
  
  for (mouseID in names(Body_motion)) {
    Bouts[[mouseID]] <- list()
    Bouts[[mouseID]][["Bouts"]] <- list()
    
    for (Bout_n in seq_along(Body_motion[[mouseID]][["Bouts"]])) {
      start_time <- Body_motion[[mouseID]][["Bouts"]][[Bout_n]]$Time[1]
      
      if (start_time > min_time && start_time < max_time) {
        Bouts[[mouseID]][["Bouts"]][[length(Bouts[[mouseID]][["Bouts"]]) + 1]] <- Body_motion[[mouseID]][["Bouts"]][[Bout_n]]
      }
    }
  }
  
  
  #Time loc these bouts few seconds before and after to their occurence
  #n_Seconds <- 2
  
  for (mouseID in names(Bouts)) {
    Bouts[[mouseID]][["Time_Range"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      Bout <- Bouts[[mouseID]][["Bouts"]][[Bout_n]]
      Start_Time <- Bout$Time[1] 
      Pre_Time <- Start_Time - n_Seconds 
      Post_Time <- Start_Time + n_Seconds
      Peak_Time <- Bout[which.max(Bout$Time),]$Time
      Peak_Pre_Time <- Peak_Time - n_Seconds
      Peak_Post_Time <- Peak_Time + n_Seconds
      
      Bouts[[mouseID]][["Time_Range"]][[Bout_n]] <- data.frame(Start_Time = Start_Time, 
                                                               Pre_Time = Pre_Time, 
                                                               Post_Time = Post_Time,
                                                               Peak_Time = Peak_Time,
                                                               Peak_Pre_Time = Peak_Pre_Time,
                                                               Peak_Post_Time = Peak_Post_Time) 
    }
  }
  
  
  #Time loc the dLight traces
  for (mouseID in names(Bouts)) {
    Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]] <- list()
    Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      Start_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Start_Time
      Pre_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Pre_Time
      Post_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Post_Time
      Peak_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Time
      Peak_Pre_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Pre_Time
      Peak_Post_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Post_Time
      
      dLight_DFF <- Body_motion[[mouseID]][["dLight_DFF"]]
      dLight_DFF <- dLight_DFF[dLight_DFF$Time > Pre_Time & dLight_DFF$Time < Post_Time,]
      dLight_DFF$Time_Loc <-  dLight_DFF$Time - Start_Time
      Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]][[Bout_n]] <- dLight_DFF
      
      dLight_DFF <- Body_motion[[mouseID]][["dLight_DFF"]]
      dLight_DFF <- dLight_DFF[dLight_DFF$Time > Peak_Pre_Time & dLight_DFF$Time < Peak_Post_Time,]
      dLight_DFF$Time_Loc <-  dLight_DFF$Time - Peak_Time
      Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]][[Bout_n]] <- dLight_DFF
      
    }
  }
  
  #Since the T times for Time loc are not exactly the same, interpolate the deltaF/F values at consensus times
  for (mouseID in names(Bouts)) {
    
    Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]] <- list()
    Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      
      Time <- seq(from = -n_Seconds, to = n_Seconds, by = 0.02)
      data <- Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]][[Bout_n]]
      Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]] <- as.data.frame(approx(x = data$Time_Loc, y = data[,curve_col_index], xout = Time))
      colnames(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]]) <- c("Time", "dLight_DFF")
      
      Time <- seq(from = -n_Seconds, to = n_Seconds, by = 0.02)
      data <- Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]
      Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]] <- as.data.frame(approx(x = data$Time_Loc, y = data[,curve_col_index], xout = Time))
      colnames(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]) <- c("Time", "dLight_DFF")
    }
  }
  
  Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]]
  
  #Combine all the interpolated traces
  All_traces <- list()
  counter <- 1
  for (mouseID in names(Bouts)) {
    for (Bout_n in seq_along(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]])) {
      All_traces[[counter]] <- Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]]
      All_traces[[counter]]$TraceID <- as.character(rep(counter, times = nrow(All_traces[[counter]])))
      All_traces[[counter]]$mouseID <- rep(mouseID, times = nrow(All_traces[[counter]]))
      counter <- counter + 1
    }
  }
  
  from_BoutStart <- bind_rows(All_traces)
  from_BoutStart <- from_BoutStart %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  All_traces <- list()
  counter <- 1 
  for (mouseID in names(Bouts)) {
    for (Bout_n in seq_along(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]])) {
      All_traces[[counter]] <- Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]
      All_traces[[counter]]$TraceID <- as.character(rep(counter, times = nrow(All_traces[[counter]])))
      All_traces[[counter]]$mouseID <- rep(mouseID, times = nrow(All_traces[[counter]]))
      counter <- counter + 1
    }
  }
  
  from_PeakSpeed <- bind_rows(All_traces)
  from_PeakSpeed <- from_PeakSpeed %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  
  plot <- from_BoutStart |>
    tidyplot(x = Time, y = dLight_DFF, color = Group) |>
    add_mean_line() |>
    add_ci95_ribbon() |>
    adjust_title(title) |>
    adjust_x_axis_title("Time from bout onset (sec.)") |>
    adjust_x_axis(limits = c(-10, 10), breaks = seq(from = -10, to = 10, by = 1)) |>
    adjust_y_axis(limits = ylim, breaks = seq(from = -10, to = 10, by = 1)) |>
    adjust_y_axis_title("dLight dF/F (%)") #|>
  
  ggsave(paste(title, "_BoutOnset.pdf", sep = ""),
         width = 4, 
         height = 4, 
         plot = plot)
  
  plot
  
  #Collect the number of bouts per mouse
  nBouts <- c()
  for (mouseID_n in 1:length(Bouts)) {
    nBouts[mouseID_n] <- length(Bouts[[mouseID_n]][["Bouts"]])
  }
  
  nBouts_df <- data.frame(mouseID = names(Bouts), nBouts = nBouts)
  nBouts_df <- nBouts_df %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  print(nBouts_df)
  write.csv(nBouts_df, paste(title, "nBouts.csv", sep = "_")) #save the data to CSV
  
  nBouts_df_pergroup <- aggregate(nBouts ~ Group, data = nBouts_df, FUN = sum, na.rm = TRUE)
  print(nBouts_df_pergroup)
  write.csv(nBouts_df_pergroup, paste(title, "nBoutsperGroup.csv", sep = "_")) #save the data to CSV
  
  
}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peak to Bout Correlation")
Plot_PeakSpeed(curve_to_timeloc = "AIN01", n_Seconds = 5, ylim = c(-10,7), min_time = 5, max_time = 300, title = "Pre-DCZ")
Plot_PeakSpeed(curve_to_timeloc = "AIN01", n_Seconds = 5, ylim = c(-10,7), min_time = 500, max_time = 1500, title = "Post-DCZ")
Plot_PeakSpeed(curve_to_timeloc = "AIN01", n_Seconds = 5, ylim = c(-10,7), min_time = 1500, max_time = 2100, title = "Pushed")
Plot_PeakSpeed(curve_to_timeloc = "DFF_LOESSsub", n_Seconds = 2, ylim = c(-6,5), min_time = 5, max_time = 300, title = "Pre-DCZ_SubstractedfromLOESS")
Plot_PeakSpeed(curve_to_timeloc = "DFF_LOESSsub", n_Seconds = 2, ylim = c(-6,5), min_time = 500, max_time = 1500, title = "Post-DCZ_SubstractedfromLOESS")
Plot_PeakSpeed(curve_to_timeloc = "DFF_LOESSsub", n_Seconds = 2, ylim = c(-6,5), min_time = 1500, max_time = 2090, title = "Pushed_SubstractedfromLOESS")
Plot_PeakSpeed(curve_to_timeloc = "DFF_LOESSsub", n_Seconds = 10, ylim = c(-6,5), min_time = 5, max_time = 300, title = "Pre-DCZ_SubstractedfromLOESS_10sec")
Plot_PeakSpeed(curve_to_timeloc = "DFF_LOESSsub", n_Seconds = 10, ylim = c(-6,5), min_time = 500, max_time = 1500, title = "Post-DCZ_SubstractedfromLOESS_10sec")
Plot_PeakSpeed(curve_to_timeloc = "DFF_LOESSsub", n_Seconds = 10, ylim = c(-6,5), min_time = 1500, max_time = 2090, title = "Pushed_SubstractedfromLOESS_10sec")

Plot_BoutStart(curve_to_timeloc = "AIN01", n_Seconds = 3, ylim = c(-10,7), min_time = 5, max_time = 300, title = "Pre-DCZ Injection")
Plot_BoutStart(curve_to_timeloc = "AIN01", n_Seconds = 3, ylim = c(-10,7), min_time = 500, max_time = 1500, title = "Post-DCZ Injection")
Plot_BoutStart(curve_to_timeloc = "AIN01", n_Seconds = 3, ylim = c(-10,7), min_time = 1500, max_time = 2100, title = "Pushed")

#------------------------------------------------------------------------
#Generate an average trace per mouse for dLight prior and post-bout occurence and bout peak-speed
Extract_Bouts_TimeLoc <- function(n_Seconds, min_time, max_time) {
  
  Bouts <- list()
  
  for (mouseID in names(Body_motion)) {
    Bouts[[mouseID]] <- list()
    Bouts[[mouseID]][["Bouts"]] <- list()
    
    for (Bout_n in seq_along(Body_motion[[mouseID]][["Bouts"]])) {
      start_time <- Body_motion[[mouseID]][["Bouts"]][[Bout_n]]$Time[1]
      
      if (start_time > min_time && start_time < max_time) {
        Bouts[[mouseID]][["Bouts"]][[length(Bouts[[mouseID]][["Bouts"]]) + 1]] <- Body_motion[[mouseID]][["Bouts"]][[Bout_n]]
      }
    }
  }
  
  
  #Time loc these bouts few seconds before and after to their occurence
  #n_Seconds <- 2
  
  for (mouseID in names(Bouts)) {
    Bouts[[mouseID]][["Time_Range"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      Bout <- Bouts[[mouseID]][["Bouts"]][[Bout_n]]
      Start_Time <- Bout$Time[1] 
      Pre_Time <- Start_Time - n_Seconds 
      Post_Time <- Start_Time + n_Seconds
      Peak_Time <- Bout[which.max(Bout$Time),]$Time
      Peak_Pre_Time <- Peak_Time - n_Seconds
      Peak_Post_Time <- Peak_Time + n_Seconds
      
      Bouts[[mouseID]][["Time_Range"]][[Bout_n]] <- data.frame(Start_Time = Start_Time, 
                                                               Pre_Time = Pre_Time, 
                                                               Post_Time = Post_Time,
                                                               Peak_Time = Peak_Time,
                                                               Peak_Pre_Time = Peak_Pre_Time,
                                                               Peak_Post_Time = Peak_Post_Time) 
    }
  }
  
  
  #Time loc the dLight traces
  for (mouseID in names(Bouts)) {
    Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]] <- list()
    Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      Start_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Start_Time
      Pre_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Pre_Time
      Post_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Post_Time
      Peak_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Time
      Peak_Pre_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Pre_Time
      Peak_Post_Time <- Bouts[[mouseID]][["Time_Range"]][[Bout_n]]$Peak_Post_Time
      
      dLight_DFF <- Body_motion[[mouseID]][["dLight_DFF"]]
      dLight_DFF <- dLight_DFF[dLight_DFF$Time > Pre_Time & dLight_DFF$Time < Post_Time,]
      dLight_DFF$Time_Loc <-  dLight_DFF$Time - Start_Time
      Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]][[Bout_n]] <- dLight_DFF
      
      dLight_DFF <- Body_motion[[mouseID]][["dLight_DFF"]]
      dLight_DFF <- dLight_DFF[dLight_DFF$Time > Peak_Pre_Time & dLight_DFF$Time < Peak_Post_Time,]
      dLight_DFF$Time_Loc <-  dLight_DFF$Time - Peak_Time
      Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]][[Bout_n]] <- dLight_DFF
      
    }
  }
  
  #Since the T times for Time loc are not exactly the same, interpolate the deltaF/F values at consensus times
  for (mouseID in names(Bouts)) {
    
    Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]] <- list()
    Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]] <- list()
    
    for (Bout_n in seq_along(Bouts[[mouseID]][["Bouts"]])) {
      
      Time <- seq(from = -n_Seconds, to = n_Seconds, by = 0.02)
      data <- Bouts[[mouseID]][["TimeLoc_dLight(BoutStart)"]][[Bout_n]]
      Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]] <- as.data.frame(approx(x = data$Time_Loc, y = data[,curve_col_index], xout = Time))
      colnames(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(BoutStart)"]][[Bout_n]]) <- c("Time", "dLight_DFF")
      
      Time <- seq(from = -n_Seconds, to = n_Seconds, by = 0.02)
      data <- Bouts[[mouseID]][["TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]
      Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]] <- as.data.frame(approx(x = data$Time_Loc, y = data[,curve_col_index], xout = Time))
      colnames(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]) <- c("Time", "dLight_DFF")
    }
  }
  
  return(Bouts)
  
}

#Generate a function for plotting average deltaF/F curve per mouse and not per bout
Plot_Average_Trace <- function(Object, Time_Loc_List, Title, x_Title) {
  
  Bouts <- Object
  
  for (mouseID in names(Bouts)) {
    
    if (length(Bouts[[mouseID]][["Bouts"]]) > 0) {
      
      #Extract all the traces and store them in a dataframe by binding the columns
      Traces <- Bouts[[mouseID]][[Time_Loc_List]] %>% 
        map(~ .x["dLight_DFF"]) %>%
        bind_cols()
      
      #Average the traces together at each T time  
      Average_dLightDFF_Trace <- rowMeans(Traces)
      Traces <- NULL
      Time <- Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[1]]$Time
      Bouts[[mouseID]][["AverageTrace_dLightDFF(PeakSpeed)"]] <- data.frame(Time = Time,
                                                                            Average_dLightDFF_Trace = Average_dLightDFF_Trace,
                                                                            mouseID = rep(mouseID, times = length(Time)))
    }
  }
  
  Average_Traces <- Bouts %>%
    map(~ .x[["AverageTrace_dLightDFF(PeakSpeed)"]]) %>%
    bind_rows(.id = "mouseID")
  
  #Add a group column
  Average_Traces <- Average_Traces %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  plot <- Average_Traces |>
    tidyplot(x = Time, y = Average_dLightDFF_Trace, color = Group) |>
    add_mean_line() |>
    add_ci95_ribbon() |>
    adjust_title(Title) |>
    adjust_x_axis_title(x_Title) |>
    #adjust_x_axis(limits = c(0.5, 5.5), breaks = seq(from = 0.5, to = 5, by = 1)) |>
    adjust_x_axis(limits = c(-1.5, 1), breaks = seq(from = -1.5, to = 1, by = 0.5)) |>
    adjust_y_axis(limits = c(-12, 12), breaks = seq(from = -10, to = 10, by = 2.5)) |>
    adjust_y_axis_title("dLight dF/F (%)")
  
  #Save the plot
  ggsave(
    filename = paste(Title, ".pdf", sep = ""), 
    plot = plot,         # Saves the last plot displayed
    #device = "png",             # Essential for transparency
    bg = "transparent",         # Removes the white background
    width = 6,                  # Width in inches
    height = 5,                 # Height in inches
    units = "in",
    dpi = 300,                  # High resolution for publication
    limitsize = TRUE
  )  
}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peak to Bout Correlation")

Plot_Average_Trace(Object = Extract_Bouts_TimeLoc(n_Seconds = 1.5,
                                                  min_time = 5,
                                                  max_time = 300),
                   Time_Loc_List = "Interpolated_TimeLoc_dLight(BoutStart)",
                   Title = "AverageTrace_PreDCZ_BoutOnset_PerMouse",
                   x_Title = "Time from bout occurence (sec.)")

Plot_Average_Trace(Object = Extract_Bouts_TimeLoc(n_Seconds = 1.5,
                                                  min_time = 500,
                                                  max_time = 1500),
                   Time_Loc_List = "Interpolated_TimeLoc_dLight(BoutStart)",
                   Title = "AverageTrace_PostDCZ_BoutOnset_PerMouse",
                   x_Title = "Time from bout occurence (sec.)")

Plot_Average_Trace(Object = Extract_Bouts_TimeLoc(n_Seconds = 1.5,
                                                  min_time = 1500,
                                                  max_time = 2100),
                   Time_Loc_List = "Interpolated_TimeLoc_dLight(BoutStart)",
                   Title = "AverageTrace_PostDCZ(ForcedLocomotion)_BoutOnset_PerMouse",
                   x_Title = "Time from bout occurence (sec.)")

Plot_Average_Trace(Object = Extract_Bouts_TimeLoc(n_Seconds = 1.5,
                                                  min_time = 5,
                                                  max_time = 300),
                   Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
                   Title = "AverageTrace_PreDCZ_BoutPeakSpeed_PerMouse",
                   x_Title = "Time from peak speed (sec.)")

Plot_Average_Trace(Object = Extract_Bouts_TimeLoc(n_Seconds = 1.5,
                                                  min_time = 500,
                                                  max_time = 1500),
                   Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
                   Title = "AverageTrace_PostDCZ_BoutPeakSpeed_PerMouse",
                   x_Title = "Time from peak speed (sec.)")

Plot_Average_Trace(Object = Extract_Bouts_TimeLoc(n_Seconds = 1.5,
                                                  min_time = 1500,
                                                  max_time = 2100),
                   Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
                   Title = "AverageTrace_PostDCZ(ForcedLocomotion)_BoutPeakSpeed_PerMouse",
                   x_Title = "Time from peak speed (sec.)")


#------------------------------------------------------------------------
#Plot density line 
Density_Line <- function(Object, Time_Loc_List, Group, nBins, Color_Option, Title, x_Title, Density_Range){
  
  Bouts <- Object
  All_traces <- list()
  counter <- 1
  for (mouseID in names(Bouts)) {
    for (Bout_n in seq_along(Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]])) {
      All_traces[[counter]] <- Bouts[[mouseID]][["Interpolated_TimeLoc_dLight(PeakSpeed)"]][[Bout_n]]
      All_traces[[counter]]$TraceID <- as.character(rep(counter, times = nrow(All_traces[[counter]])))
      All_traces[[counter]]$mouseID <- rep(mouseID, times = nrow(All_traces[[counter]]))
      counter <- counter + 1
    }
  }
  
  from_BoutStart <- bind_rows(All_traces)
  from_BoutStart <- from_BoutStart %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  All_traces <- list()
  counter <- 1 
  for (mouseID in names(Bouts)) {
    for (Bout_n in seq_along(Bouts[[mouseID]][[Time_Loc_List]])) {
      All_traces[[counter]] <- Bouts[[mouseID]][[Time_Loc_List]][[Bout_n]]
      All_traces[[counter]]$TraceID <- as.character(rep(counter, times = nrow(All_traces[[counter]])))
      All_traces[[counter]]$mouseID <- rep(mouseID, times = nrow(All_traces[[counter]]))
      counter <- counter + 1
    }
  }
  
  All_traces <- bind_rows(All_traces)
  All_traces <- All_traces %>%
    mutate(
      Group = case_when(
        mouseID  == "J0777" ~ "fDIO-mCherry",
        mouseID  == "J0782" ~ "fDIO-mCherry",
        mouseID  == "J0859" ~ "fDIO-mCherry",
        mouseID  == "K0871" ~ "fDIO-mCherry",
        
        mouseID  == "J0779" ~ "CoffFon-hM4Di",
        mouseID  == "J0784" ~ "CoffFon-hM4Di",
        mouseID  == "K0331" ~ "CoffFon-hM4Di",
        mouseID  == "K0994" ~ "CoffFon-hM4Di",
        mouseID  == "K0861" ~ "CoffFon-hM4Di",
        
        mouseID  == "K0029" ~ "ConFon-hM4Di",
        mouseID  == "K0863" ~ "ConFon-hM4Di",
        mouseID  == "K0865" ~ "ConFon-hM4Di",
        mouseID  == "K0990" ~ "ConFon-hM4Di"
      )
    )
  
  #Rename the traceID
  All_traces$TraceID <- paste(All_traces$mouseID,All_traces$TraceID, sep = "_")
  
  #Subset the group to plot
  All_traces <- All_traces[All_traces$Group == Group,]
  
  #Plot
plot <- ggplot(data = All_traces, aes(x = Time, y = dLight_DFF, group = TraceID)) + 
    stat_line_density(aes(color = after_stat(density), fill = after_stat(density)), 
                      bins = nBins, drop = FALSE, na.rm = TRUE) +
    #scale_color_viridis_c(option = "magma") +
    scale_fill_viridis_c(option = Color_Option, limits = Density_Range, values = c(0,0.5,1)) +
    scale_x_continuous(limits = c(-2,2), expand = c(0, 0), breaks = seq(from = -2, to = 2, by = 0.5)) + 
    scale_y_continuous(limits = c(-25,25), expand = c(0, 0), breaks = seq(from = -20, to = 20 , by = 5)) +
    theme_classic() +
    labs(
      title = Title, # Main title
      x = x_Title, # X-axis label
      y = "dLight deltaF/F (%)", # Y-axis label
    ) +
    theme(plot.title = element_text(size = 12, hjust = 0.5))

ggsave(
  filename = paste(Title, "_DensityLine.png", sep = ""), 
  plot = plot,         # Saves the last plot displayed
  #device = "png",             # Essential for transparency
  bg = "transparent",         # Removes the white background
  width = 6,                  # Width in inches
  height = 5,                 # Height in inches
  units = "in",
  dpi = 1200,                  # High resolution for publication
  limitsize = TRUE)
}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peak to Bout Correlation")
Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 5,
                                            max_time = 300), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "fDIO-mCherry",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PreDCZ_fDIO-mCherry",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 500,
                                            max_time = 1500), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "fDIO-mCherry",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PostDCZ_fDIO-mCherry",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 1500,
                                            max_time = 2100), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "fDIO-mCherry",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PostDCZ(Forced)_fDIO-mCherry",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 5,
                                            max_time = 300), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "CoffFon-hM4Di",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PreDCZ_CoffFon-hM4Di",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 500,
                                            max_time = 1500), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "CoffFon-hM4Di",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PostDCZ_CoffFon-hM4Di",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 1500,
                                            max_time = 2100), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "CoffFon-hM4Di",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PostDCZ(Forced)_CoffFon-hM4Di",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 5,
                                            max_time = 300), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "ConFon-hM4Di",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PreDCZ_ConFon-hM4Di",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 500,
                                            max_time = 1500), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "ConFon-hM4Di",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PostDCZ_ConFon-hM4Di",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

Density_Line(Object = Extract_Bouts_TimeLoc(n_Seconds = 2,
                                            min_time = 1500,
                                            max_time = 2100), 
             Time_Loc_List = "Interpolated_TimeLoc_dLight(PeakSpeed)",
             Group = "ConFon-hM4Di",
             nBins = 100,
             Color_Option = "turbo",
             Title = "PostDCZ(Forced)_ConFon-hM4Di",
             x_Title = "Time from bout peak speed (sec.)",
             Density_Range = NULL)

#------------------------------------------------------------------------
#Correlate the velocity and dLight DFF at every T times
velocity_dLightDFF_datapoints <- list()

for (mouseID in names(Body_motion)) {
  #Extract the velocity and dLight_DFF data for each frame
  velocity_df <- Body_motion[[mouseID]][["Sleap_data"]][,c("Time", "Torso.velocity")]
  dLightDFF_df <- Body_motion[[mouseID]][["dLight_DFF"]]
  
  #Approximate velocity on each T time of the dLight data
  approximation <- approx(x = velocity_df$Time, y = velocity_df$Torso.velocity, xout = dLightDFF_df$Time)
  dLightDFF_df$Torso.velocity <- approximation$y
  
  #Add a column for the mouseID
  dLightDFF_df$mouseID <- rep(mouseID, times = nrow(dLightDFF_df))
  
  #Add the dataframe to the list
  velocity_dLightDFF_datapoints[[mouseID]] <- dLightDFF_df
  
  #Clean the objects generated
  velocity_df <- NULL
  dLightDFF_df <- NULL
  approximation <- NULL
}

#Stack the dataframes
velocity_dLightDFF_datapoints <- dplyr::bind_rows(velocity_dLightDFF_datapoints)

#Add a group column
velocity_dLightDFF_datapoints <- velocity_dLightDFF_datapoints %>%
  mutate(
    Group = case_when(
      mouseID  == "J0777" ~ "fDIO-mCherry",
      mouseID  == "J0782" ~ "fDIO-mCherry",
      mouseID  == "J0859" ~ "fDIO-mCherry",
      mouseID  == "K0871" ~ "fDIO-mCherry",
      
      mouseID  == "J0779" ~ "CoffFon-hM4Di",
      mouseID  == "J0784" ~ "CoffFon-hM4Di",
      mouseID  == "K0331" ~ "CoffFon-hM4Di",
      mouseID  == "K0994" ~ "CoffFon-hM4Di",
      mouseID  == "K0861" ~ "CoffFon-hM4Di",
      
      mouseID  == "K0029" ~ "ConFon-hM4Di",
      mouseID  == "K0863" ~ "ConFon-hM4Di",
      mouseID  == "K0865" ~ "ConFon-hM4Di",
      mouseID  == "K0990" ~ "ConFon-hM4Di"
    )
  )


Velocity_DFF_DensityPlot <- function(Start_time, 
                                     End_time, 
                                     Group, 
                                     Title, 
                                     density_limits
                                     ) {
  Toplot <- velocity_dLightDFF_datapoints[velocity_dLightDFF_datapoints$Time > Start_time &
                                            velocity_dLightDFF_datapoints$Time < End_time &
                                            velocity_dLightDFF_datapoints$Group == Group,]
  
plot <-  ggplot(Toplot, mapping = aes(x = Torso.velocity, y = AIN01)) +
    #eom_point(alpha = 1/10) +
    geom_pointdensity(size = 0.1, 
                      adjust = 5) +
    scale_color_viridis(option = "magma", 
                        limits = density_limits,
                        oob = scales::squish, #Prevents the color from disappearing if saturated
                        alpha = 1) + 
    xlim(0,1000) +
    ylim(-8,12) +
    theme_classic(base_size = 10) + 
    #geom_hline(yintercept = 0) + 
    #geom_vline(xintercept = 0) +
    labs(
      title = Title, # Main title
      x = "Velocity (mm/sec.)", # X-axis label
      y = "dLight deltaF/F (%)", # Y-axis label
    ) +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) + #Centre the title
    annotate("text", x = 800, y = 12, label = paste(as.character(nrow(Toplot)), "Frames"), color = "gray40", size = 3) #Add the number of frames

#Save the plot
ggsave(
  filename = paste(Title, "_Velocity&dLightDFF_DensityPlot.png", sep = ""), 
  plot = plot,         # Saves the last plot displayed
  #device = "png",             # Essential for transparency
  bg = "transparent",         # Removes the white background
  width = 6,                  # Width in inches
  height = 5,                 # Height in inches
  units = "in",
  dpi = 300,                  # High resolution for publication
  limitsize = TRUE
)  

}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peak to Bout Correlation")
Velocity_DFF_DensityPlot(Start_time = 5, 
                         End_time = 300, 
                         Group = "fDIO-mCherry", 
                         Title = "fDIO-mCherry (PreDCZ)", 
                         density_limits = c(0, 0.005)
                         )

Velocity_DFF_DensityPlot(Start_time = 500, 
                         End_time = 1500, 
                         Group = "fDIO-mCherry", 
                         Title = "fDIO-mCherry (PostDCZ)", 
                         density_limits = c(0, 0.005)
                         )

Velocity_DFF_DensityPlot(Start_time = 1500, 
                         End_time = 2100, 
                         Group = "fDIO-mCherry", 
                         Title = "fDIO-mCherry (PostDCZ + Forced Locomotion)", 
                         density_limits = c(0, 0.005)
                         )

Velocity_DFF_DensityPlot(Start_time = 5, 
                         End_time = 300, 
                         Group = "CoffFon-hM4Di", 
                         Title = "CoffFon-hM4Di (PreDCZ)", 
                         density_limits = c(0, 0.005)
                         )

Velocity_DFF_DensityPlot(Start_time = 500, 
                         End_time = 1500, 
                         Group = "CoffFon-hM4Di", 
                         Title = "CoffFon-hM4Di (PostDCZ)", 
                         density_limits = c(0, 0.005)
                         )

Velocity_DFF_DensityPlot(Start_time = 1500, 
                         End_time = 2100, 
                         Group = "CoffFon-hM4Di", 
                         Title = "CoffFon-hM4Di (PostDCZ + Forced Locomotion)", 
                         density_limits = c(0, 0.005)
                         )

Velocity_DFF_DensityPlot(Start_time = 5, 
                         End_time = 300, 
                         Group = "ConFon-hM4Di", 
                         Title = "ConFon-hM4Di (PreDCZ)", 
                         density_limits = c(0, 0.005)
)

Velocity_DFF_DensityPlot(Start_time = 500, 
                         End_time = 1500, 
                         Group = "ConFon-hM4Di", 
                         Title = "ConFon-hM4Di (PostDCZ)", 
                         density_limits = c(0, 0.005)
)

Velocity_DFF_DensityPlot(Start_time = 1500, 
                         End_time = 2100, 
                         Group = "ConFon-hM4Di", 
                         Title = "ConFon-hM4Di (PostDCZ + Forced Locomotion)", 
                         density_limits = c(0, 0.005)
)

#------------------------------------------------------------------------
#Plot sum deltaF/F and velocity per 10 sec., labelling PreDCZ, Transition, PostDCZ, PostDCZ+Forced
#Interpolate dLight deltaDFF 
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["Sleap_data"]]$dLightDFF <- approx(Body_motion[[mouseID]][["dLight_DFF"]]$Time, 
                                                             Body_motion[[mouseID]][["dLight_DFF"]]$AIN01,
                                                             Body_motion[[mouseID]][["Sleap_data"]]$Time)$y
  Body_motion[[mouseID]][["Velocity/dLight"]] <- data.frame(Time = Body_motion[[mouseID]][["Sleap_data"]]$Time,
                                                            Velocity = Body_motion[[mouseID]][["Sleap_data"]]$Torso.velocity,
                                                            dLightDFF = Body_motion[[mouseID]][["Sleap_data"]]$dLightDFF)
}

#Assign the time values to PreDCZ, Transition, PostDCZ, PostDCZ_Forced
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["Velocity/dLight"]] <- Body_motion[[mouseID]][["Velocity/dLight"]] %>%
    mutate(
      Period = case_when(
        Time  >= 0 & Time < 300 ~ "PreDCZ",
        Time  >= 300 & Time < 500 ~ "Transition",
        Time  >= 500 & Time < 1500 ~ "PostDCZ",
        Time  >= 1500 & Time < 2100 ~ "PostDCZ (Forced Locomotion)"
      )
    )
}


#Assign the bins
for (mouseID in names(Body_motion)) {
  bin_time <- 1 #in sec.
  nbins <- 2100 / bin_time
  
  Body_motion[[mouseID]][["Velocity/dLight"]]$BinID <-   cut(Body_motion[[mouseID]][["Velocity/dLight"]]$Time, 
                                                             breaks = seq(from = 0, to = 2100, by = bin_time),
                                                             labels = 1:nbins)
  
}

#Average the dLight and velocity values inside the bins
for (mouseID in names(Body_motion)) {
  temp1 <- Body_motion[[mouseID]][["Velocity/dLight"]] %>%
    group_by(BinID) %>%
    summarize(mean_velocity = mean(Velocity))
  
  temp2 <- Body_motion[[mouseID]][["Velocity/dLight"]] %>%
    group_by(BinID) %>%
    summarize(SD_dLightDFF = sd(dLightDFF))

Body_motion[[mouseID]][["Mean_Velocity/dLight"]] <- data.frame(BinID = temp1$BinID,
                                                                 mean_velocity = temp1$mean_velocity,
                                                                 SD_dLightDFF = temp2$SD_dLightDFF)

#Reassign the period time
Body_motion[[mouseID]][["Mean_Velocity/dLight"]]$BinID <- as.numeric(Body_motion[[mouseID]][["Mean_Velocity/dLight"]]$BinID)
Body_motion[[mouseID]][["Mean_Velocity/dLight"]] <- Body_motion[[mouseID]][["Mean_Velocity/dLight"]] %>%
  mutate(
    Period = case_when(
      BinID  >= 0 & BinID < 300/bin_time ~ "PreDCZ",
      BinID  >= 300/bin_time & BinID < 500/bin_time ~ "Transition",
      BinID  >= 500/bin_time & BinID < 1500/bin_time ~ "PostDCZ",
      BinID  >= 1500/bin_time & BinID <= 2100/bin_time ~ "PostDCZ (Forced Locomotion)")
  )

#Add the mouseID
Body_motion[[mouseID]][["Mean_Velocity/dLight"]]$mouseID <- rep(mouseID, times = nrow(Body_motion[[mouseID]][["Mean_Velocity/dLight"]]))

}

#Stack the dataframes together
Binned_velocity_dLightDFF <- Body_motion %>%
  map(~ .x[["Mean_Velocity/dLight"]]) %>%
  bind_rows()

#Add group column
Binned_velocity_dLightDFF <- Binned_velocity_dLightDFF %>%
  mutate(
    Group = case_when(
      mouseID  == "J0777" ~ "fDIO-mCherry",
      mouseID  == "J0782" ~ "fDIO-mCherry",
      mouseID  == "J0859" ~ "fDIO-mCherry",
      mouseID  == "K0871" ~ "fDIO-mCherry",
      
      mouseID  == "J0779" ~ "CoffFon-hM4Di",
      mouseID  == "J0784" ~ "CoffFon-hM4Di",
      mouseID  == "K0331" ~ "CoffFon-hM4Di",
      mouseID  == "K0994" ~ "CoffFon-hM4Di",
      mouseID  == "K0861" ~ "CoffFon-hM4Di",
      
      mouseID  == "K0029" ~ "ConFon-hM4Di",
      mouseID  == "K0863" ~ "ConFon-hM4Di",
      mouseID  == "K0865" ~ "ConFon-hM4Di",
      mouseID  == "K0990" ~ "ConFon-hM4Di"
    )
  )

dLightSDvsMeanVelocity_PerTime <- function(Group) {
  Toplot <- Binned_velocity_dLightDFF[Binned_velocity_dLightDFF$Group == Group,]
  Toplot <- Toplot[sample(nrow(Toplot)), ] #Randomize the rows appearance
  
  #Remove the NAs
  Toplot <- Toplot[!is.na(Toplot$Period),]
  Toplot$Time <- Toplot$BinID*bin_time
  ggplot(Toplot, mapping = aes(x = mean_velocity, y = SD_dLightDFF, color = Time)) +
    theme_classic() +
    xlim(0,300) +
    ylim(0,6) +
    labs(
      title = Group, # Main title
      x = "Velocity Average (mm/s per 5s)", # X-axis label
      y = "dLight deltaF/F Standard Deviation (per 5s)", # Y-axis label
    ) +
    theme(plot.title = element_text(size = 12, hjust = 0.5)) + #Centre the title
    geom_point(alpha = 0.9, size = 1) + 
    scale_color_viridis_c(option = "turbo", values = c(0, 0.2, 0.3, 0.7, 1), breaks = c(0, 300, 500, 1500, 2100))

}

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Fibrephotometry/dLight_DCZinj5min/Output/Peak to Bout Correlation") 
dLightSDvsMeanVelocity_PerTime(Group = "fDIO-mCherry")
dLightSDvsMeanVelocity_PerTime(Group = "CoffFon-hM4Di")
dLightSDvsMeanVelocity_PerTime(Group = "ConFon-hM4Di")




dLightSDvsMeanVelocity_PerPeriod <- function(Group) {
  Toplot <- Binned_velocity_dLightDFF[Binned_velocity_dLightDFF$Group == Group,]
  Toplot <- Toplot[sample(nrow(Toplot)), ] #Randomize the rows appearance
  
  #Remove the NAs
  Toplot <- Toplot[!is.na(Toplot$Period),]
  Toplot$Time <- Toplot$BinID*bin_time
  
  #Transform the Period to factor so they appear in the temporal order
  Toplot$Period <- factor(Toplot$Period, levels = c("PreDCZ", "Transition", "PostDCZ", "PostDCZ (Forced Locomotion)"))
  color <- viridis(100, option = "turbo")
  
  #Discard the transition period
  Toplot <- Toplot %>%
    dplyr::filter(Period %in% c("PreDCZ", "PostDCZ", "PostDCZ (Forced Locomotion)"))
  
  
  p <- ggplot(Toplot, mapping = aes(x = mean_velocity, y = SD_dLightDFF, color = Period)) +
    theme_classic() + 
    stat_ellipse(type = "norm", linetype = "dashed", linewidth = 1.5, level = 0.95) +
    xlim(-100,250) +
    ylim(-1,6) +
    labs(
      title = Group, # Main title
      x = "Velocity Average (mm/s)", # X-axis label
      y = "dLight deltaF/F Standard Deviation", # Y-axis label
    ) +
    theme(plot.title = element_text(size = 12, hjust = 0.5), #Centre the title
          legend.position = "bottom",         # Move legend to bottom
          legend.direction = "vertical"     # Make it lay out side-by-side        
    ) + 
    geom_point(alpha = 0.9, size = 0.5) + 
    #scale_color_viridis_d(option = "cividis")
    scale_color_manual(values = c("PreDCZ" = color[5], 
                                  "PostDCZ" = color[64],
                                  "PostDCZ (Forced Locomotion)" = color[98])) +
    theme(legend.position = "none")
  
  ggMarginal(p, type = "density", groupFill = TRUE, alpha = 0.99)
  
}

dLightSDvsMeanVelocity_PerPeriod(Group = "fDIO-mCherry")
dLightSDvsMeanVelocity_PerPeriod(Group = "CoffFon-hM4Di")
dLightSDvsMeanVelocity_PerPeriod(Group = "ConFon-hM4Di")


#------------------------------------------------------------------------
#Plot average velocity and dLightDFF across the phases
#First remove the rows with NA
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["Velocity/dLight"]] <- na.omit(Body_motion[[mouseID]][["Velocity/dLight"]])
}

#Then compute the average velocity per period
for (mouseID in names(Body_motion)) {
  temp1 <- Body_motion[[mouseID]][["Velocity/dLight"]] %>%
    group_by(Period) %>%
    summarize(mean_velocity = mean(Velocity))
  
  temp2 <- Body_motion[[mouseID]][["Velocity/dLight"]] %>%
    group_by(Period) %>%
    summarize(SD_dLightDFF = sd(dLightDFF))
  
  Body_motion[[mouseID]][["Mean_Velocity/dLight"]] <- data.frame(Period = temp1$Period,
                                                                 mean_velocity = temp1$mean_velocity,
                                                                 SD_dLightDFF = temp2$SD_dLightDFF)
}

#Add the mouseID
for (mouseID in names(Body_motion)) {
  Body_motion[[mouseID]][["Mean_Velocity/dLight"]]$mouseID <- rep(mouseID, times = nrow(Body_motion[[mouseID]][["Mean_Velocity/dLight"]]))
}


#Stack the dataframes together
Mean_velocity_dLightDFF <- Body_motion %>%
  map(~ .x[["Mean_Velocity/dLight"]]) %>%
  bind_rows()

#Add group column
Mean_velocity_dLightDFF <- Mean_velocity_dLightDFF %>%
  mutate(
    Group = case_when(
      mouseID  == "J0777" ~ "fDIO-mCherry",
      mouseID  == "J0782" ~ "fDIO-mCherry",
      mouseID  == "J0859" ~ "fDIO-mCherry",
      mouseID  == "K0871" ~ "fDIO-mCherry",
      
      mouseID  == "J0779" ~ "CoffFon-hM4Di",
      mouseID  == "J0784" ~ "CoffFon-hM4Di",
      mouseID  == "K0331" ~ "CoffFon-hM4Di",
      mouseID  == "K0994" ~ "CoffFon-hM4Di",
      mouseID  == "K0861" ~ "CoffFon-hM4Di",
      
      mouseID  == "K0029" ~ "ConFon-hM4Di",
      mouseID  == "K0863" ~ "ConFon-hM4Di",
      mouseID  == "K0865" ~ "ConFon-hM4Di",
      mouseID  == "K0990" ~ "ConFon-hM4Di"
    )
  )

#Plot
library(ggpubr)
Mean_velocity_dLightDFF$Period <-  factor(Mean_velocity_dLightDFF$Period, levels = c("PreDCZ", "Transition", "PostDCZ", "PostDCZ (Forced Locomotion)"))
Mean_velocity_dLightDFF$Group <-  factor(Mean_velocity_dLightDFF$Group, levels = c("fDIO-mCherry", "CoffFon-hM4Di", "ConFon-hM4Di"))

#Remove the transition time
Mean_velocity_dLightDFF <- Mean_velocity_dLightDFF %>%
  dplyr::filter(!Period %in% "Transition")

#Specify the sex
Mean_velocity_dLightDFF <- Mean_velocity_dLightDFF %>%
  mutate(Sex = case_when(
    mouseID == "J0782" ~ "M",
    mouseID == "J0859" ~ "M",
    mouseID == "J0777" ~ "F",
    mouseID == "J0784" ~ "M",
    mouseID == "J0779" ~ "F",
    mouseID == "J0782" ~ "M",
    mouseID == "K0029" ~ "F",
    mouseID == "K0992" ~ "M",
    mouseID == "K0871" ~ "F",
    mouseID == "K0331" ~ "M",
    mouseID == "K0994" ~ "M",
    mouseID == "K0863" ~ "M",
    mouseID == "K0865" ~ "M",
    mouseID == "J0868" ~ "F",
    mouseID == "K0990" ~ "F",
    mouseID == "K0861" ~ "F"
    
  ))


Mean_velocity_dLightDFF |>
  tidyplot(x = Period, y = mean_velocity, color = Group) |>
  add_mean_bar(alpha = 0.3) |>
  add_sem_errorbar() |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  add_test_asterisks(hide_info = TRUE) |>
  adjust_y_axis(title = "Mean Velocity (mm/s)")


Mean_velocity_dLightDFF |>
  tidyplot(x = Period, y = SD_dLightDFF, color = Group) |>
  add_mean_bar(alpha = 0.3) |>
  add_sem_errorbar() |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  add_test_asterisks(hide_info = TRUE) |>
  adjust_y_axis(title = "dLight deltaF/F Standard Deviation")

#Run two-way ANOVA 

#####Velocity##########

#Check assumption
#Normality
ggqqplot(Mean_velocity_dLightDFF, "mean_velocity", ggtheme = theme_bw()) +
  facet_grid(Group ~ Period)
#Note: All the points fall approximately along the reference line, for each cell. So we can assume normality of the data.

#Homogeneity of variance assumption
Mean_velocity_dLightDFF %>% levene_test(mean_velocity ~ Group*Period)
#The Levene’s test not significant so we can assume homogeneity

#ANOVA
res.aov <- Mean_velocity_dLightDFF %>% anova_test(mean_velocity ~ Group*Period)
res.aov
#ANOVA has a significant interaction, so I will proceed to Tukey pairwise comparison
pwc <- Mean_velocity_dLightDFF %>% 
  group_by(Period) %>%
  tukey_hsd(mean_velocity ~ Group) 

write.csv(pwc, "MeanVelocity_PerPeriod_twowayANOVATukeyHSD.csv")

#####deltaF/F##########

#Check assumption
#Normality
ggqqplot(Mean_velocity_dLightDFF, "SD_dLightDFF", ggtheme = theme_bw()) +
  facet_grid(Group ~ Period)
#Note: All the points fall approximately along the reference line, for each cell. So we can assume normality of the data.

#Homogeneity of variance assumption
Mean_velocity_dLightDFF %>% levene_test(SD_dLightDFF ~ Group*Period)
#The Levene’s test not significant so we can assume homogeneity

#ANOVA
res.aov <- Mean_velocity_dLightDFF %>% anova_test(SD_dLightDFF ~ Group*Period)
res.aov
#ANOVA has a significant interaction, so I will proceed to Tukey pairwise comparison
pwc <- Mean_velocity_dLightDFF %>% 
  group_by(Period) %>%
  tukey_hsd(SD_dLightDFF ~ Group) 

write.csv(pwc, "SDdeltaFF_PerPeriod_twowayANOVATukeyHSD.csv")

