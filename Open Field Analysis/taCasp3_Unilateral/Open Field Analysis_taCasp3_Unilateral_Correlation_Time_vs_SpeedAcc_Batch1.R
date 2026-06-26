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

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1")
#------------------------------------------------------------------------------
#Get the image resolution in cm

#Load the image
image_filenames <- list.files(pattern="*.png")
img <- readPNG(image_filenames[1])
h <- nrow(img)
w <- ncol(img)
windows(width = 10, height = 10 * (h/w)) 
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, type = "n", xlim = c(0, w), ylim = c(h, 0), asp = 1, axes = FALSE)
rasterImage(img, 0, h, w, 0)
pts <- locator(2) #Add two points on the corner of the openfield

#Show the points
points(pts$x, pts$y, col = "red", pch = 19, cex = 1.5)
x1 <- pts$x[1]
y1 <- pts$y[1]
x2 <- pts$x[2]
y2 <- pts$y[2]
Pixel_Distance <- sqrt((x2-x1)^2+(y2-y1)^2)

Resolution <- 39.5/Pixel_Distance
Resolution

#Note: The value I got for this batch one pixel is equivalent to 0.0448cm, or 0.448mm
Resolution <- 0.448

#-----------------------------------------------------------------------------

#Read the CSV files
#List the csv files within the document
filenames <- list.files(pattern="*.csv")
nfiles <- length(filenames)

filenames
#-----------------------------------------------------------------------------
#First generate the velocity dataframe to store all the data
Velocity_Data <- data.frame(matrix(nrow = nfiles, ncol = 16)) 
colnames(Velocity_Data) <- c("Group", 
                             "MouseID",
                             "Mode Velocity (mm/s)",
                             "% Time in Mobility",
                             "Number of Bouts", 
                             "Average Bout Length (sec)",
                             "Average Peak Velocity (mm/s)",
                             "Average Peak Acceleration (mm/s2)",
                             "Average Peak Deceleration (mm/s2)",
                             "Average Peak Delay (sec)",
                             "Average Percent Time to Reach Peak Speed",
                             "Total Distance Travelled_0-3min (mm)",
                             "Total Distance Travelled_3-6min (mm)",
                             "Total Distance Travelled_6-9min (mm)",
                             "Total Distance Travelled_9-12min (mm)",
                             "Total Distance Travelled_12-15min (mm)"
)
#Generate a list to store all the bouts for each animals
Bouts_Amplitude <- vector(mode='list', length=nfiles)
names(Bouts_Amplitude) <- filenames

for (n in 1:nfiles) {
  
  print(paste("Processing", filenames[n]))
  
  setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1")
  data <- as.data.frame(read.csv(filenames[n]))
  data <- t(data)
  
  #Set another working directory to store all the output files
  setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1/Output")
  
  #-----------------------------------------------------------------------------
  #Extract the information about the ID of the mouse
  Group <- sapply(strsplit(filenames, "_"), `[`, 1)
  MouseID <- sapply(strsplit(filenames, "_"), `[`, 3)
  MouseID <- gsub('.csv','',MouseID)
  
  #Store the group and ID to the velocity dataframe
  Velocity_Data[n,"Group"] <- Group[n]
  Velocity_Data[n,"MouseID"] <- MouseID[n]
  #-----------------------------------------------------------------------------
  #Reannotate the first rows of the dataframe
  temp <- paste(data[-1,1], "_", data[-1,2], sep = "")
  data <- data[,3:ncol(data)]
  data <- data[-1,]
  
  data <- cbind(temp, data)
  data <- t(data)
  colnames(data) <- data[1,]
  data <- data[-1,]
  rownames(data) <- c(1:nrow(data))
  
  #Convert each frame by its actual time value in seconds
  #Note: The camera has a number of 29.83 frames per second, so one frame is:
  frame_time <- 1/29.83
  
  frame <- c(0)
  for (x in 2:nrow(data)) {
    frame[x] <- x * frame_time
  }
  
  
  
  data <- as.data.frame(data) #Convert everything to dataframe
  data <- data.frame(apply(data, 2, function(x) as.numeric(as.character(x)))) #Make sure all values are numeric
  
  #------------------------------------------------------------------------
  #Subset the body data and filter the data, making sure that each X and Y positions have a likelihood score >0.9
  back <- as.data.frame(cbind(frame, data[13:15]))
  back<- back[back$back_likelihood>0.9,]
  
  #Convert the x y coordinates into mm scale
  back[,2:3] <- back[,2:3] * Resolution
  
  #------------------------------------------------------------------------
  #Downsample the data to the time resolution
  time_resolution <- 0.1 #Time resolution in second
  
  Ref_frames <- seq(from = 0, to = 900, by = time_resolution)
  Closest_frames <- Closest(back$frame, Ref_frames, which = FALSE, na.rm = FALSE)
  
  New_frames <- c()
  for (x in 1:length(Closest_frames)) {
    New_frames[x] <- Closest_frames[[x]]
  }
  
  New_frames <- unique(New_frames) #Some frames have been assigned in double, so need to take only the unique ones
  
  Frame_indices <- c()
  for (x in 1:length(New_frames)) {
    Frame_indices[x] <- which(back$frame == New_frames[x])
  }
  
  back <- back[Frame_indices,]
  
  #------------------------------------------------------------------------
  #Extract the values for velocity
  velocity <- function(X1, X2, Y1, Y2, T1, T2) { 
    temp <- (sqrt((X2 - X1)^2 + (Y2-Y1)^2))/(T2-T1)
    return (temp)
  }
  
  speed_back <- c(0)
  for (x in 2:nrow(back)) {
    speed_back[x] <- velocity(X1 = back$back_x[x-1],
                              X2 = back$back_x[x],
                              Y1 = back$back_y[x-1],
                              Y2 = back$back_y[x],
                              T1 = back$frame[x-1],
                              T2 = back$frame[x])
    
  }
  
  #------------------------------------------------------------------------
  #Extract the values for acceleration
  acceleration <- function(X1, X2, T1, T2) { 
    temp <- (X2 - X1)/(T2-T1)
    return (temp)
  }
  acceleration_back <- c(0)
  for (x in 2:length(speed_back)) {
    acceleration_back[x] <- acceleration(X1 = speed_back[x-1],
                                         X2 = speed_back[x],
                                         T1 = back$frame[x-1],
                                         T2 = back$frame[x])
    
  }
  
  # create data
  data <- data.frame(back$frame,speed_back, acceleration_back)
  colnames(data) <- c("Frame", "Speed", "Acceleration")
  data <- data[data$Frame < 900, ]
  
  #------------------------------------------------------------------------
  # Extract the bouts 
  #Note: Bouts will be defined as speed is defined as an occurence where the speed is suddenly increased >25mm/s for at least 600ms 
  
  #Set the treshold corresponding to 25mm/s
  Treshold <- 25
  
  #Extract the bouts
  temp <- data$Speed
  temp[temp < Treshold] <- 0 #Assign to 0 all speed values lower to 0 to extract the bouts
  data$OriginalSpeed <- data$Speed #Keep an original speed column
  data$Speed <- temp
  
  index <- which(data$Speed %in% c(0))
  
  Start_Index <- c()
  End_Index <- c()
  
  for (x in 1:(length(index)-1)) {
    x1 <- index[x] #Get the index of the start of the potential bout
    x2 <- index[x+1] #Get the index of the end of the potential bout
    
    temp <- data[x1:x2,] #Extract the dataframe at that bout
    
    Start <- temp[1, "Frame"]
    End_600ms <- Start + 0.6 #Estimate where the 600ms ends
    
    End <- Closest(temp$Frame, End_600ms, which = FALSE, na.rm = FALSE) #Get the closest value to that end of 600ms
    
    End_index <- which(temp$Frame %in% End) #Get the index of that closest value
    
    temp <- temp[2:End_index,] #First row will have a speed of 0, extract all the other rows within that 600ms
    
    if (all(temp$Speed > 25) == TRUE) {
      Start_Index[x] <- x1
      End_Index[x] <- x2
    } 
  }
  
  Bouts <- cbind(Start_Index, End_Index)
  Bouts <- na.omit(Bouts) #Remove the NA to extract only the bouts indexes
  Bouts <- as.data.frame(Bouts)

  #------------------------------------------------------------------------
  #Extract the time of bout occurence and the max velocity, acceleration and deceleration
  Time <- c()
  MaxSpeed <- c()
  MaxAcceleration <- c()
  MaxDeceleration <- c()
  
  for (x in 1:nrow(Bouts)) {
    Time[x] <- data[Bouts$Start_Index[x],]$Frame
    MaxSpeed[x] <- max(data[Bouts$Start_Index[x]:Bouts$End_Index[x],]$Speed)
    MaxAcceleration[x] <- max(data[Bouts$Start_Index[x]:Bouts$End_Index[x],]$Acceleration)
    MaxDeceleration[x] <- min(data[Bouts$Start_Index[x]:Bouts$End_Index[x],]$Acceleration)
  }

  Bouts_Amplitude[[n]] <- data.frame(Time, MaxSpeed, MaxAcceleration, MaxDeceleration)
}

saveRDS(Bouts_Amplitude, "BoutsAmplitude_Batch1.rds")

#------------------------------------------------------------------------
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1/Output")

#Read the lists and combine them together
Batch1 <- readRDS("BoutsAmplitude_Batch1.rds")
Batch2 <- readRDS("BoutsAmplitude_Batch2.rds")
Batch3 <- readRDS("BoutsAmplitude_Batch3.rds")

Bouts_Amplitude <- c(Batch1, Batch2, Batch3)

#Add a Group column
Animals <- names(Bouts_Amplitude)

for (x in 1:length(Bouts_Amplitude)) {
  Bouts_Amplitude[[x]]$Group <- rep(str_split(Animals, '_')[[x]][1], times = nrow(Bouts_Amplitude[[x]]))
}
#------------------------------------------------------------------------
Bouts_Correlation <- bind_rows(Bouts_Amplitude)

Bouts_Correlation$Group <- factor(Bouts_Correlation$Group, levels = c("NotaCasp3", "CoffFontaCasp3"))


ggplot(Bouts_Correlation, aes(Time, MaxSpeed, colour = Group)) +
  geom_point(shape = 16, size = 2, alpha = 0.25) +
  geom_smooth(method = "lm", se = FALSE, size =1) +
  stat_cor(method = "pearson", 
           aes(label = paste(..rr.label.., sep = "~`,`~")),
           label.x = 700, 
           label.y = c(540, 500, 460), 
           size = 4, 
           show.legend = FALSE) +
  scale_color_manual(values = c("NotaCasp3" = "#000000", "CoffFontaCasp3" = "#FF0700")) +
  labs(x = "Time (s)", y = "Peak Velocity (mm/s)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 900)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 600)) +
  theme_classic() +
  theme(legend.position = "none")

ggplot(Bouts_Correlation, aes(Time, MaxAcceleration, colour = Group)) +
  geom_point(shape = 16, size = 2, alpha = 0.25) +
  geom_smooth(method = "lm", se = FALSE, size =1) +
  stat_cor(method = "pearson", 
           aes(label = paste(..rr.label.., sep = "~`,`~")),
           label.x = 700, 
           label.y = c(2400, 2200, 2000), 
           size = 4, 
           show.legend = FALSE) +
  scale_color_manual(values = c("NotaCasp3" = "#000000", "CoffFontaCasp3" = "#FF0700")) +
  labs(x = "Time (s)", y = "Peak Acceleration (mm/s2)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 900)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) +
  theme_classic() +
  theme(legend.position = "none")

ggplot(Bouts_Correlation, aes(Time, MaxDeceleration, colour = Group)) +
  geom_point(shape = 16, size = 2, alpha = 0.25) +
  geom_smooth(method = "lm", se = FALSE, size =1) +
  stat_cor(method = "pearson", 
           aes(label = paste(..rr.label.., sep = "~`,`~")),
           label.x = 700, 
           label.y = c(-2400, -2200, -2000), 
           size = 4, 
           show.legend = FALSE) +
  scale_color_manual(values = c("NotaCasp3" = "#000000", "CoffFontaCasp3" = "#FF0700")) +
  labs(x = "Time (s)", y = "Peak Deceleration (mm/s2)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 900)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, -2500)) +
  theme_classic() +
  theme(legend.position = "none")

#Calculate the pearson correlation
Corr_to_time <- as.data.frame(matrix(nrow = 3, ncol = 3))
colnames(Corr_to_time) <- unique(Bouts_Correlation$Group)
rownames(Corr_to_time) <- colnames(Bouts_Correlation[2:4])  
pValue <- Corr_to_time

for (x in 1:ncol(Corr_to_time)) {
  temp <- Bouts_Correlation[Bouts_Correlation$Group == colnames(Corr_to_time)[x],]
  
  Corr_to_time["MaxSpeed", x] <- cor(temp$Time, temp$MaxSpeed)
  Corr_to_time["MaxAcceleration", x] <- cor(temp$Time, temp$MaxAcceleration)
  Corr_to_time["MaxDeceleration", x] <- cor(temp$Time, temp$MaxDeceleration)
}


for (x in 1:ncol(pValue)) {
  temp <- Bouts_Correlation[Bouts_Correlation$Group == colnames(pValue)[x],]
  
  pValue["MaxSpeed", x] <- cor.test(temp$Time, temp$MaxSpeed)[["p.value"]]
  pValue["MaxAcceleration", x] <- cor.test(temp$Time, temp$MaxAcceleration)[["p.value"]]
  pValue["MaxDeceleration", x] <- cor.test(temp$Time, temp$MaxDeceleration)[["p.value"]]
}

write.csv(Corr_to_time, "PearsonCoeff_Time_vs_VelAcc_taCasp3Uni")
write.csv(pValue, "pValueCorr_Time_vs_VelAcc_taCasp3Uni")
Corr_to_time
pValue
