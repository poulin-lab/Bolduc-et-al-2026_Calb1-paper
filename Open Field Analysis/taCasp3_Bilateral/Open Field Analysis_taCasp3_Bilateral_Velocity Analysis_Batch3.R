library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(DescTools)
library(magick)
library(imager)
library(lattice)
library(viridis)
library(paletteer)
library(gplots)
library(sommer)

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch3")
#------------------------------------------------------------------------------
#Get the image resolution in cm

#Load the image
#image_filenames <- list.files(pattern="*.png")
#image_filenames
#img_cimg <- load.image(image_filenames[1])  # imager format
#plot(img_cimg)  # Show image

#Place two dots from the left and right side corners of the chamber, which is 39.5cm
#pts <- locator(2)  # User clicks two points

#Extract the distance in pixels
#x1 <- pts$x[1]
#y1 <- pts$y[1]
#x2 <- pts$x[2]
#y2 <- pts$y[2]

#Pixel_Distance <- sqrt((x2-x1)^2+(y2-y1)^2)
#Resolution <- 39.5/Pixel_Distance
#Resolution

#Note: The value I got for this batch one pixel is equivalent to 0.0380cm, or 0.380mm
Resolution <- 0.380

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
All_Bouts <- vector(mode='list', length=nfiles)

for (n in 1:nfiles) {

  setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch3")
  data <- as.data.frame(read.csv(filenames[n]))
data <- t(data)

#Set another working directory to store all the output files
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch3/Output")

#-----------------------------------------------------------------------------
#Extract the information about the ID of the mouse
Group <- sapply(strsplit(filenames, "_"), `[`, 1)
MouseID <- sapply(strsplit(filenames, "_"), `[`, 2)
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
#Note: Bouts will be defined as speed is defined as an occurence where the speed is suddenly increased at 25mm/s and >25mm/s for at least 400ms 
#Note: The bout starts at the moment it transition from 2x the mode, and ends when it comes back to that value

#Extract the x and y values from the kernel density graph
g <- plot(data %>%   
            ggplot( aes(x=Speed)) +
            geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
            labs(y= "Density", x = "Velocity (mm/s)") +
            ylim(0,0.21) +
            xlim(0,500))

p <- ggplot_build(g)

#Extract the mode speed    
Mode_Index <- which.max(p$data[[1]]$density)
Mode_Index <- as.numeric(Mode_Index)
Density_Speed <- p$data[[1]]$x
Mode <- Density_Speed[Mode_Index]

#Store the mode velocity to the velocity_data dataframe
Velocity_Data[n,"Mode Velocity (mm/s)"] <- Mode

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
  End_400ms <- Start + 0.4 #Estimate where the 400ms ends
  
  End <- Closest(temp$Frame, End_400ms, which = FALSE, na.rm = FALSE) #Get the closest value to that end of 200ms
  
  End_index <- which(temp$Frame %in% End) #Get the index of that closest value
  
  temp <- temp[2:End_index,] #First row will have a speed of 0, extract all the other rows within that 200ms
  
  if (all(temp$Speed > Treshold) == TRUE) {
    Start_Index[x] <- x1
    End_Index[x] <- x2
  } 
}

Bouts <- cbind(Start_Index, End_Index)
Bouts <- na.omit(Bouts) #Remove the NA to extract only the bouts indexes
Bouts <- as.data.frame(Bouts)

#------------------------------------------------------------------------
#View the bouts
#Bout_number <- 1
#Bout <- data[Bouts[Bout_number,1]:Bouts[Bout_number,2],]
#Bout

#ggplot(Bout, aes(x=Frame, y=Speed)) +
#labs(y= "Velocity (mm/s)", x = "Time (s)") +
#ylim(0,500) +
#geom_line()

#------------------------------------------------------------------------
#Add the bouts to the data object, assigning to 0 all velocities and acceleration not within the bouts 

Bouts_Velocity <- data

Bouts_Velocity[1:Bouts[1, "Start_Index"], "OriginalSpeed"] <- -0.1
Bouts_Velocity[Bouts[nrow(Bouts), "End_Index"]:nrow(Bouts_Velocity), "OriginalSpeed"] <- -0.1

for (x in 1:(nrow(Bouts)-1)) {
  
  Bouts_Velocity[Bouts[x, "End_Index"]:Bouts[x+1, "Start_Index"], "OriginalSpeed"] <- -0.1

}

New_data <- data

New_data$Group <- rep("", nrow(data))
Bouts_Velocity$Group <- rep("Bout", nrow(Bouts_Velocity))
New_data <- rbind(New_data, Bouts_Velocity)

#------------------------------------------------------------------------
# Plot the graphs for velocity and acceleration courses and density
pdf(file = paste("Velocity", "_", Group[n], "_", MouseID[n], ".pdf", sep=""),
    width = 10, # The width of the plot in inches
    height = 2) # The height of the plot in inches

plot(ggplot(New_data, aes(x=Frame, y=OriginalSpeed, group=Group, colour=Group)) +
  scale_color_manual(values=c("black","#CC6666")) +
  labs(y= "Velocity (mm/s)", x = "Time (s)") +
  ylim(0,500) +
  geom_line(size=0.01))

dev.off()


pdf(file = paste("Acceleration", "_", Group[n], "_",MouseID[n], ".pdf", sep=""),
    width = 10, # The width of the plot in inches
    height = 2) # The height of the plot in inches

plot(ggplot(data, aes(x=Frame, y=Acceleration)) +
  labs(y= "Acceleration (mm/s2)", x = "Time (s)") +
  ylim(-3000,3000) +
  geom_line(size=0.01))

dev.off()


# Make the histogram
pdf(file = paste("Acceleration Density", "_", Group[n], "_",MouseID[n], ".pdf", sep=""),
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

plot(data %>%   
  ggplot( aes(x=Acceleration)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  labs(y= "Density", x = "Acceleration (mm/s2)") +
  ylim(0,0.018) +
  xlim(-2500,2500))

dev.off()

pdf(file = paste("Velocity Density", "_", Group[n], "_",MouseID[n], ".pdf", sep=""),
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

plot(data %>%   
  ggplot( aes(x=Speed)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  labs(y= "Density", x = "Velocity (mm/s)") +
  ylim(0,0.21) +
  xlim(0,500))

dev.off()

#------------------------------------------------------------------------
#%Time immobile

#Get all the timelenght of each bouts
Time_mobile <- 0
for (x in 1:(nrow(Bouts))) {
X1 <- Bouts[x,"Start_Index"]
X2 <- Bouts[x, "End_Index"]

Time_mobile <- Time_mobile + (data[X2, "Frame"] - data[X1, "Frame"])
  }

#Extract the percentage of time mobile
Percent_Time_Mobile <- (Time_mobile/data[nrow(data),"Frame"]) *100

#Store the data in the dataframe
Velocity_Data[n,"% Time in Mobility"] <- Percent_Time_Mobile
#------------------------------------------------------------------------
#Number of bouts
Number_bouts <- nrow(Bouts)

#Store the number of bouts to the dataframe
Velocity_Data[n, "Number of Bouts"] <- Number_bouts
#------------------------------------------------------------------------
#Bouts length
Bouts_length <- c()
for (x in 1:nrow(Bouts)) {
  Start <- Bouts[x,"Start_Index"]
  End <- Bouts[x,"End_Index"]
  Bouts_length[x] <- data[End,"Frame"] - data[Start,"Frame"]
}

#Store the bouts lenght in the data.frame
Velocity_Data[n, "Average Bout Length (sec)"] <- mean(Bouts_length)

#------------------------------------------------------------------------
#Peak Velocity
Bouts_amplitude <- c()

for (x in 1:nrow(Bouts)) {
  Start <- Bouts[x,"Start_Index"]
  End <- Bouts[x,"End_Index"]
  Bouts_amplitude[x] <- max(data[Start:End,"Speed"])
}

Velocity_Data[n, "Average Peak Velocity (mm/s)"] <- mean(Bouts_amplitude)

#------------------------------------------------------------------------
#Peak Acceleration
Peak_Acceleration <- c()

for (x in 1:nrow(Bouts)) {
  Start <- Bouts[x,"Start_Index"]
  End <- Bouts[x,"End_Index"]
  Peak_Acceleration[x] <- max(data[Start:End,"Acceleration"])
}

Velocity_Data[n, "Average Peak Acceleration (mm/s2)"] <- mean(Peak_Acceleration)

#------------------------------------------------------------------------
#Peak Deceleration
Peak_Deceleration <- c()

for (x in 1:nrow(Bouts)) {
  Start <- Bouts[x,"Start_Index"]
  End <- Bouts[x,"End_Index"]
  Peak_Deceleration[x] <- min(data[Start:End,"Acceleration"])
}

Velocity_Data[n, "Average Peak Deceleration (mm/s2)"] <- mean(Peak_Deceleration)

#------------------------------------------------------------------------
#Time to reach peak velocity
Peak_delay <- c()
Peak_Index <- c()

for (x in 1:nrow(Bouts)) {
  Start <- Bouts[x,"Start_Index"]
  End <- Bouts[x,"End_Index"]
  Peak_Index[x] <- Start + (which.max(data[Start:End,"Speed"])-1)
  Peak_delay[x] <- data[Peak_Index[x], "Frame"] - data[Start, "Frame"]
}

#Store the data within the velocity dataframe
Velocity_Data[n, "Average Peak Delay (sec)"] <- mean(Peak_delay)

#------------------------------------------------------------------------
#Extract the velocities values one second before and three seconds after the movement initiation
Bouts_list <- vector(mode='list', length=Number_bouts)

#Note: For that analysis, we can't get the speed from the bouts before the first second, and after the three last seconds
min_index <- nrow(data[data$Frame <1,])
max_index <- nrow(data) - nrow(data[data$Frame >897,])

#Filter the bouts based on these min and max index found up there
Filtered_Bouts <- Bouts[Bouts$Start_Index > min_index & Bouts$End_Index < max_index,]

for (x in 1:nrow(Filtered_Bouts)) {
  
  Initiation_Index <- Filtered_Bouts[x,"Start_Index"]
  Start <- data[Initiation_Index, "Frame"]
  
  MinusOne <- Start - 1
  PlusThree <- Start + 3
  
  MinusOne <- Closest(data$Frame, MinusOne, which = FALSE, na.rm = FALSE)
  PlusThree <- Closest(data$Frame, PlusThree, which = FALSE, na.rm = FALSE)
  
  MinusOne_Index <- which(data$Frame == MinusOne)
  PlusThree_Index <- which(data$Frame == PlusThree)
  
  #Extract the speed and acceleration values from -1 and +3 seconds after the peak acceleration
  temp <- data[MinusOne_Index:PlusThree_Index, ]

  new_temp <- temp
  
  for  (y in 1:nrow(temp)) {
  #Recenter all the frame values to the start of the bout
  new_temp[y,"Frame"] <- temp[y,"Frame"] - temp[((Initiation_Index-MinusOne_Index)+1),"Frame"]
  
  
  }

  Bouts_list[[x]] <- new_temp
  
}

All_Bouts[[n]] <- Bouts_list #Store the bouts within the list
names(All_Bouts)[n] <- paste(Group[n], MouseID[n], sep = "_") #Rename the list element with the corresponding mouse group and ID

#Discard the NULL bouts generated. I can't find why they are generated
All_Bouts[[n]] <- All_Bouts[[n]] %>% discard(is.null)
#------------------------------------------------------------------------
#%Time to reach peak speed within the bouts
Perc_Delay_Peak <- (Peak_delay/Bouts_length)*100

#Store the data within the velocity dataframe
Velocity_Data[n, "Average Percent Time to Reach Peak Speed"] <- mean(Perc_Delay_Peak)

#------------------------------------------------------------------------
#Total Distance Travelled in 3 minutes timelapses
First_3min <- back[back$frame < 180,]
Second_3min <- back[back$frame > 180 & back$frame < 360,]
Third_3min <- back[back$frame > 360 & back$frame < 540,]
Fourth_3min <- back[back$frame > 540 & back$frame < 720,]
Fifth_3min <- back[back$frame > 720 & back$frame < 900,]

#Define a function to measure the total distance travelled within the timeplases
Distance_Travelled <- function(Timelapse) { 
  temp <- 0
  for (x in 1:(nrow(Timelapse)-1)) {
    temp <- temp + sqrt((Timelapse$back_x[x+1] - Timelapse$back_x[x])^2 + (Timelapse$back_y[x+1]-Timelapse$back_y[x])^2)
  }
  return(temp)
}

Velocity_Data[n, "Total Distance Travelled_0-3min (mm)"] <- Distance_Travelled(Timelapse = First_3min)
Velocity_Data[n, "Total Distance Travelled_3-6min (mm)"] <- Distance_Travelled(Timelapse = Second_3min)
Velocity_Data[n, "Total Distance Travelled_6-9min (mm)"] <- Distance_Travelled(Timelapse = Third_3min)
Velocity_Data[n, "Total Distance Travelled_9-12min (mm)"] <- Distance_Travelled(Timelapse = Fourth_3min)
Velocity_Data[n, "Total Distance Travelled_12-15min (mm)"] <- Distance_Travelled(Timelapse = Fifth_3min)



}#It's the bracket for line 59

#Export the results to CSV
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

write.csv(Velocity_Data,"Velocity_Data_Batch3.csv", row.names = FALSE)

#------------------------------------------------------------------------
#Save all the bouts
saveRDS(All_Bouts, file = "All_Bouts_DREADD_Batch3.rds")

View(All_Bouts)

#------------------------------------------------------------------------

#Extract the all the bouts for all the groups

#Generate a function to extract all the bouts for each group
Group_All_Bouts <- function(Group) { 
  #Extract the groups index in the All_Bouts list
  Bouts_Group_Index <- which(startsWith(names(All_Bouts), Group), arr.ind = FALSE, useNames = TRUE)
  
  #Extract the bouts of the appropriate length
  Combined_Bouts <- c()
  for  (x in Bouts_Group_Index) {
    for  (y in 1:length(All_Bouts[[x]])) {
      temp <- All_Bouts[[x]][[y]]
      
      #Extract only the bouts where there are no missing frames, which is 41 frames at the time resolution we have
      if (nrow(temp) == 41) {
        temp <- temp$OriginalSpeed
        Combined_Bouts <- rbind(Combined_Bouts, temp)
      } 
    }
  }
  
  Combined_Bouts <- as.data.frame(Combined_Bouts)
  return(Combined_Bouts)  
  
}

#Then, extract the all the bouts for all the groups and store them within the Groups_Bouts list
Groups <- c("ConFon-taCasp3", "NotaCasp3")
Groups_Bouts <- vector(mode='list', length(Groups))
names(Groups_Bouts) <- Groups

for  (x in 1:length(Groups)) {
  Groups_Bouts[[x]] <- Group_All_Bouts(Group = Groups[x])
}

#------------------------------------------------------------------------
#Plot the average bout curves

#Generate a data frame to store the average bout per group
Average_Bout_Curve <- data.frame(matrix(nrow = 41, ncol = 4)) 
colnames(Average_Bout_Curve) <- c("Frame", Groups)
Average_Bout_Curve$Frame <- All_Bouts[[1]][[1]]$Frame

for  (x in 1:length(Groups)) {
  temp <- Groups_Bouts[[x]]
  Average_Bout_Curve[,x+1] <- colMeans(temp)
}

Average_Bout_Curve <- data.frame(rep(Average_Bout_Curve[,1], times=3), 
                                 c(Average_Bout_Curve[,2],
                                   Average_Bout_Curve[,3],
                                   Average_Bout_Curve[,4]),
                                 c(
                                   rep("ConFon-taCasp3", times = 41),
                                   rep("NotaCasp3", times = 41)))

colnames(Average_Bout_Curve) <- c("Frame", "Velocity_Value", "Group")


plot(ggplot(Average_Bout_Curve, aes(x=Frame, y=Velocity_Value, color = Group)) +
       labs(y= "Velocity (mm/s)", x = "Time (s)") +
       ylim(0,80) +
       geom_line(size=1.5))



#------------------------------------------------------------------------
#Generate a heatmap of all the bouts velocity by creating a function
Bouts_HeatMap <- function(Group, Random_Bouts_Number) { 
  Combined_Bouts <- Groups_Bouts[[Group]]
  
  #Randomly sample some bouts
  Combined_Bouts <- Combined_Bouts[sample(nrow(Combined_Bouts), Random_Bouts_Number), ] 
  
  #Order the randomly picked bouts by their peak amplitude 1 second after the initiation of it
  max <- as.vector(apply(Combined_Bouts[,10:20], MARGIN = 1, FUN = max))
  max_index <- sort(max, decreasing =  TRUE, index.return = TRUE)
  
  New_Combined_Bouts <- c()
  for  (x in 1:nrow(Combined_Bouts)) {
    temp <- Combined_Bouts[max_index$ix[x],]
    New_Combined_Bouts <- rbind(New_Combined_Bouts, temp)
    
  } 
  
  New_Combined_Bouts <- t(as.matrix(New_Combined_Bouts))
  
  #Plot a heatmap of the bouts -1s and +3s after the bout initiation
  levelplot(New_Combined_Bouts, col.regions = paletteer_c("viridis::inferno", 100),
            cuts = 100, pretty = TRUE, labels = TRUE, at=seq(0,200,length=100), useRaster = TRUE)
  
}

Bouts_HeatMap(Group = "NotaCasp3",
              Random_Bouts_Number = 300)

Bouts_HeatMap(Group = "ConFon-taCasp3",
              Random_Bouts_Number = 300)
