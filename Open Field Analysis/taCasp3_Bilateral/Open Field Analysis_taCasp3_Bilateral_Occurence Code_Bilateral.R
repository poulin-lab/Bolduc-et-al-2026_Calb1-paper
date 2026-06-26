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

#------------------------------------------------------------------------
#Save all the bouts
#saveRDS(All_Bouts, file = "All_Bouts_DREADD_Batch1.rds")

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch1/Output")
l1 <- readRDS(file = "All_Bouts_taCasp3_Batch1.rds")
l2 <- readRDS(file = "All_Bouts_taCasp3_Batch2.rds")

#Merge the lists together
All_Bouts <- append(l1, l2)
#------------------------------------------------------------------------
#Extract the average acceleration curve
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
        temp <- temp$Acceleration
        Combined_Bouts <- rbind(Combined_Bouts, temp)
      } 
    }
  }
  
  Combined_Bouts <- as.data.frame(Combined_Bouts)
  return(Combined_Bouts)  
  
}

#Then, extract the all the bouts for all the groups and store them within the Groups_Bouts list
Groups <- c("CoffFontaCasp3", "ConFontaCasp3", "NotaCasp3")
Groups_Bouts <- vector(mode='list', length(Groups))
names(Groups_Bouts) <- Groups

for  (x in 1:length(Groups)) {
  Groups_Bouts[[x]] <- Group_All_Bouts(Group = Groups[x])
}

#Plot the bout curves
#Generate a function for the average bout curve
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
                                 c(rep("CoffFontaCasp3", times = 41),
                                   rep("ConFontaCasp3", times = 41),
                                   rep("NotaCasp3", times = 41)))

colnames(Average_Bout_Curve) <- c("Frame", "Acceleration_Value", "Group")

pdf(file = "Average Acceleration Curve.pdf",
    width = 7, # The width of the plot in inches
    height = 4) # The height of the plot in inches

ggplot(Average_Bout_Curve, aes(x=Frame, y=Acceleration_Value, color = Group)) +
  labs(y= "Acceleration (mm/s2)", x = "Time (s)") +
  ylim(-220,420) +
  scale_color_manual(values = c("CoffFontaCasp3" = "red", "ConFontaCasp3" = "chartreuse4", "NotaCasp3" = "black")) +
  geom_line(size=0.7) + theme_classic()

dev.off()

#------------------------------------------------------------------------
#Extract the average velocity curve
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
Groups <- c("CoffFontaCasp3", "ConFontaCasp3", "NotaCasp3")
Groups_Bouts <- vector(mode='list', length(Groups))
names(Groups_Bouts) <- Groups

for  (x in 1:length(Groups)) {
  Groups_Bouts[[x]] <- Group_All_Bouts(Group = Groups[x])
}

#Plot the bout curves
#Generate a function for the average bout curve
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
                                 c(rep("CoffFontaCasp3", times = 41),
                                   rep("ConFontaCasp3", times = 41),
                                   rep("NotaCasp3", times = 41)))

colnames(Average_Bout_Curve) <- c("Frame", "Velocity_Value", "Group")

pdf(file = "Average Velocity Curve.pdf",
    width = 7, # The width of the plot in inches
    height = 4) # The height of the plot in inches

ggplot(Average_Bout_Curve, aes(x=Frame, y=Velocity_Value, color = Group)) +
  labs(y= "Velocity (mm/s)", x = "Time (s)") +
  ylim(0,120) +
  scale_color_manual(values = c("CoffFontaCasp3" = "red", "ConFontaCasp3" = "chartreuse4", "NotaCasp3" = "black")) +
  geom_line(size=0.7) + theme_classic()


dev.off()

#------------------------------------------------------------------------
#Generate a heatmap of all the bouts velocity by creating a function
Bouts_HeatMap <- function(Group, Random_Bouts_Number) { 
  Combined_Bouts <- Groups_Bouts[[Group]]
  
  #Randomly sample some bouts
  Combined_Bouts <- Combined_Bouts[sample(nrow(Combined_Bouts), Random_Bouts_Number), ] 
  
  #Order the randomly picked bouts by their peak amplitude 1 second after the initiation of it
  max <- as.vector(apply(Combined_Bouts[,10:40], MARGIN = 1, FUN = mean))
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

Bouts_HeatMap(Group = "NotaCasp3", Random_Bouts_Number = 150)

Bouts_HeatMap(Group = "CoffFontaCasp3", Random_Bouts_Number = 150)

Bouts_HeatMap(Group = "ConFontaCasp3", Random_Bouts_Number = 150)

#------------------------------------------------------------------------
#Extract the frequencies of bout peak speed, duration time, acceleration and deceleration peak

Bouts_Frequency <- vector(mode='list', length(All_Bouts))

for (i in 1:length(Bouts_Frequency)) {
  
  Peak_Speed <- c()
  Duration <- c()
  Peak_Acceleration <- c()
  Peak_Deceleration <- c()
  
  for  (x in 1:length(All_Bouts[[i]])) {
    
    Bout <- All_Bouts[[i]][[x]]
    
    Bout <- Bout[Bout$Frame > 0,]
    
    #Extract the bouts that only finish within the 3 seconds following its initiation
    if (any(Bout$Speed == 0)) {
      End_Bout_Index <- which(Bout[Bout$Frame > 0,"Speed"] == 0)
      End_Bout_Index <- End_Bout_Index[1]
      Bout <- Bout[Bout$Frame > 0,]
      
      Bout_Speed <- Bout[1:End_Bout_Index, "Speed"]
      Bout_Acceleration <- Bout[1:End_Bout_Index, "Acceleration"]
      
      #Extract all the bout peak speed
      Peak_Speed[x] <- max(Bout_Speed)
      
      #Extract all the bout duration
      Duration[x] <- Bout[End_Bout_Index, "Frame"] - Bout[1, "Frame"]
      
      #Extract all the bout peak acceleration
      Peak_Acceleration[x] <- max(Bout_Acceleration)
      
      #Extract all the bout peak acceleration
      Peak_Deceleration[x] <- min(Bout_Acceleration)
      
    }
    
    #Remove the NAs generated within the Peak_Speed Matrix
    Peak_Speed <- Peak_Speed[!is.na(Peak_Speed)]
    Peak_Speed[Peak_Speed > 300] <- 300
    Duration <- Duration[!is.na(Duration)]
    Duration[Duration > 3] <- 3
    Peak_Acceleration <- Peak_Acceleration[!is.na(Peak_Acceleration)]
    Peak_Acceleration[Peak_Acceleration > 1500] <- 1500 #Assign the acceleration values higher than 2000 to 3000
    Peak_Deceleration <- Peak_Deceleration[!is.na(Peak_Deceleration)]
    Peak_Deceleration[Peak_Deceleration < -1500] <- -1500 #Assign the acceleration values lower than -2000 to 3000
    
    
  }
  
  names(Bouts_Frequency)[i] <- names(All_Bouts[i]) #Rename the list element with the corresponding mouse group and ID
  
  Peak_Speed_Hist <- hist(Peak_Speed, breaks = seq(0, 300, length.out = 31), ylim = c(0,100))
  Peak_Speed_Frequency <- data.frame(Peak_Speed_Hist$breaks[-1],(Peak_Speed_Hist$counts/sum(Peak_Speed_Hist$counts)))
  Bouts_Frequency[[i]][[1]] <- Peak_Speed_Frequency
  names(Bouts_Frequency[[i]])[[1]] <- "Peak Speed"
  
  Duration_Hist <- hist(Duration, breaks = seq(0, 3, length.out = 31), ylim = c(0,300))
  Duration_Frequency <- data.frame(Duration_Hist$breaks[-1],Duration_Hist$counts/sum(Duration_Hist$counts))
  Bouts_Frequency[[i]][[2]] <- Duration_Frequency
  names(Bouts_Frequency[[i]])[[2]] <- "Duration"
  
  Peak_Acceleration_Hist <- hist(Peak_Acceleration, breaks = seq(0, 1500, length.out = 16), ylim = c(0,100))
  Peak_Acceleration_Frequency <- data.frame(Peak_Acceleration_Hist$breaks[-1],Peak_Acceleration_Hist$counts/sum(Peak_Acceleration_Hist$counts))
  Bouts_Frequency[[i]][[3]] <- Peak_Acceleration_Frequency
  names(Bouts_Frequency[[i]])[[3]] <- "Peak Acceleration"
  
  Peak_Deceleration_Hist <- hist(Peak_Deceleration, breaks = seq(0, -1500, length.out = 16), ylim = c(0,100))
  Peak_Deceleration_Frequency <- data.frame(Peak_Deceleration_Hist$breaks[-length(Peak_Deceleration_Hist$breaks)],Peak_Deceleration_Hist$counts/sum(Peak_Deceleration_Hist$counts))
  Bouts_Frequency[[i]][[4]] <- Peak_Deceleration_Frequency
  names(Bouts_Frequency[[i]])[[4]] <- "Peak Deceleration"
  
}

Peak_Speed
#------------------------------------------------------------------------
#Generate the graph for peak speed

#First, extract the all the peak speed frequencies per group
Groups <- c("CoffFontaCasp3", "ConFontaCasp3", "NotaCasp3")
Groups_Bouts_Frequency <- vector(mode='list', length(Groups))
names(Groups_Bouts_Frequency) <- Groups

Average_Occurence <- function(Parameter, Range) {
  
  for (x in 1:length(Groups)) {
    
    #Get the group index within the Bouts_Frequency list
    Group_Index <- which(startsWith(names(Bouts_Frequency), Groups[x]), arr.ind = FALSE, useNames = TRUE)
    
    #Create an empty dataframe to store the data
    df <- data.frame(matrix(nrow = length(Group_Index), ncol = length(Bouts_Frequency[[1]][[Parameter]][,1]))) 
    colnames(df) <- Bouts_Frequency[[1]][[Parameter]][,1]
    
    for (y in Group_Index) {
      df[which(Group_Index == y),] <- Bouts_Frequency[[y]][[Parameter]][,2]
      Groups_Bouts_Frequency[[x]] <- df
    }
  }
  
  
  #Then average all the bouts frequency per group
  df <- data.frame(matrix(nrow = length(Bouts_Frequency[[1]][[Parameter]][,1])))
  df[,1] <- Bouts_Frequency[[1]][[Parameter]][,1]
  colnames(df) <- "Breaks"
  
  for (x in 1:length(Groups_Bouts_Frequency)) {
    df[,x+1]<- colMeans(Groups_Bouts_Frequency[[x]])
  }
  
  Occurence_Curve <- data.frame(rep(df[,1], times=3), 
                                c(df[,2],
                                  df[,3],
                                  df[,4]),
                                c(rep("CoffFontaCasp3", times = nrow(df)),
                                  rep("ConFontaCasp3", times = nrow(df)),
                                  rep("NotaCasp3", times = nrow(df))))
  
  colnames(Occurence_Curve) <- c("Parameter", "Occurence", "Group")
  
  
  plot(ggplot(Occurence_Curve, aes(x=Parameter, y=Occurence, color = Group)) +
         labs(y= "Occurence", x = Parameter) +
         ylim(Range) +
         geom_line(size=1.5))
  
  return(Groups_Bouts_Frequency)
  
}

#Then, export the occurence into CSV files

MouseIDs <- c("E765",
              "E772",
              "E774",
              "E767",
              "E770",
              "E801",
              "F273",
              "F275",
              "F283",
              "F284",
              "F290",
              "D696",
              "E701",
              "E768",
              "F191",
              "F272",
              "F285",
              "F286",
              "F305"
)


#Compute the cumulative occurence for Peak Speed
Peak_Speed <- Average_Occurence(Parameter = "Peak Speed", Range = c(0,1))
temp <- t(rbind(Peak_Speed[[1]],Peak_Speed[[2]],Peak_Speed[[3]]))
colnames(temp) <- MouseIDs
data <- temp
for (x in 1:nrow(data)) {
  for (y in 1:ncol(data)) {
    data[x,y] <- sum(temp[1:x,y])
  }
}
write.csv(data,"Peak_Speed_Cumulative_Occurrence.csv")

#Compute the cumulative occurence for Bouts Duration
Duration <- Average_Occurence(Parameter = "Duration", Range = c(0,1))
temp <- t(rbind(Duration[[1]],Duration[[2]],Duration[[3]]))
colnames(temp) <- MouseIDs
data <- temp
for (x in 1:nrow(data)) {
  for (y in 1:ncol(data)) {
    data[x,y] <- sum(temp[1:x,y])
  }
}
write.csv(data,"Duration_Cumulative_Occurrence.csv")

#Compute the cumulative occurence for Peak Acceleration
Peak_Acceleration <- Average_Occurence(Parameter = "Peak Acceleration", Range = c(0,0.3))
temp <- t(rbind(Peak_Acceleration[[1]],Peak_Acceleration[[2]],Peak_Acceleration[[3]]))
colnames(temp) <- MouseIDs
data <- temp
for (x in 1:nrow(data)) {
  for (y in 1:ncol(data)) {
    data[x,y] <- sum(temp[1:x,y])
  }
}
write.csv(data,"Peak_Acceleration_Cumulative_Occurrence.csv")

#Compute the cumulative occurence for Peak Deceleration
Peak_Deceleration <- Average_Occurence(Parameter = "Peak Deceleration", Range = c(0,0.3))
temp <- t(rbind(Peak_Deceleration[[1]],Peak_Deceleration[[2]],Peak_Deceleration[[3]]))
colnames(temp) <- MouseIDs
temp <- temp[rev(rownames(temp)), ] #Reverse the rows, so the cumulative occurence can start from 0 to -1500, not the opposite
data <- temp
for (x in 1:nrow(data)) {
  for (y in 1:ncol(data)) {
    data[x,y] <- sum(temp[1:x,y])
  }
}

write.csv(data,"Peak_Deceleration__Cumulative_Occurrence.csv")
#-------------------------------------------------------------------------------
library(ggplot2)

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Bilateral/Batch1/Output")

SEM <- function(x, MARGIN) {
  apply(x, MARGIN = MARGIN, FUN = sd) / sqrt(nrow(x))
}

#-------------------------------------------------------------------------------
#Generate a function to save the cumulative occurence curve and the Two-ANOVA Tukey Post-Hoc Stats

Cumulative_Occurence_Plot <- function(File, #Name for the csv file. Should end with ".csv"
                                      Plotname, #Name of the pdf file to save. Should end with ".pdf"
                                      x_title, #Character for the name of the x axis
                                      y_title, #Character for the name of the y axis
                                      xlim, #xlim for the limits in x. Should be written that way: xlim(0,1500) 
                                      ylim, #ylim for the limits in y. Should be written that way: ylim(0,1500)
                                      Group1, #Character of the group name
                                      Group1_Col_Index, #Define the column index of the group1. Should be written as; 1:4
                                      Group2, #Character of the group name
                                      Group2_Col_Index, #Define the column index of the group2. Should be written as; 5:8
                                      Group3, #Character of the group name
                                      Group3_Col_Index #Define the column index of the group2. Should be written as; 9:12
         ) 

  {

#------------------------
#Processing of the objects
  
# Load first Group
data1 <- read.csv(file = File)
X <- data1$X
data1 <- cbind(X, data1[, Group1_Col_Index])
data1$Average <- rowMeans(data1[, !names(data1) %in% "X"])

MouseIDs1 <- colnames(data1[, !names(data1) %in% c("X", "Average")])
sem1 <- SEM(data1[, MouseIDs1], MARGIN = 1)
data1$sem_min <- data1$Average - sem1
data1$sem_max <- data1$Average + sem1


# Load second Group
data2 <- read.csv(file = File)
X <- data2$X
data2 <- cbind(X, data2[, Group2_Col_Index])  
data2$Average <- rowMeans(data2[, !names(data2) %in% "X"])

MouseIDs2 <- colnames(data2[, !names(data2) %in% c("X", "Average")])
sem2 <- SEM(data2[, MouseIDs2], MARGIN = 1)
data2$sem_min <- data2$Average - sem2
data2$sem_max <- data2$Average + sem2


# Load third Group
data3 <- read.csv(file = File)
X <- data3$X
data3 <- cbind(X, data3[, Group3_Col_Index])
data3$Average <- rowMeans(data3[, !names(data3) %in% "X"])

MouseIDs3 <- colnames(data3[, !names(data3) %in% c("X", "Average")])
sem3 <- SEM(data3[, MouseIDs3], MARGIN = 1)
data3$sem_min <- data3$Average - sem3
data3$sem_max <- data3$Average + sem3

# Plot
p <- ggplot() +
  # Group1 labelled with Green ribbon 
  geom_ribbon(data = data1, aes(x = X, ymin = sem_min, ymax = sem_max), fill = "red", alpha = 0.3) +
  geom_line(data = data1, aes(x = X, y = Average), color = adjustcolor("red", alpha.f = 1), size = 0.5)
# Add individual mouse lines
#for (x in 1:length(MouseIDs1)) {
#  temp <- data.frame(X = data1$X, Average = data1[, MouseIDs1[x]])
#  p <- p + geom_line(data = temp, aes(x = X, y = Average), color = adjustcolor("darkgreen", alpha.f = 0.1), size = 1)
#}

# Group2 labelled with Red ribbon
p <- p +
  geom_ribbon(data = data2, aes(x = X, ymin = sem_min, ymax = sem_max), fill = "chartreuse4", alpha = 0.3) +
  geom_line(data = data2, aes(x = X, y = Average), color = adjustcolor("chartreuse4", alpha.f = 1), size = 0.5) 

# Add individual mouse lines
#for (x in 1:length(MouseIDs2)) {
#  temp <- data.frame(X = data2$X, Average = data2[, MouseIDs2[x]])
#  p <- p + geom_line(data = temp, aes(x = X, y = Average), color = adjustcolor("red", alpha.f = 0.1), size = 1)
#}

# Group3 labelled with Black ribbon
p <- p +
  geom_ribbon(data = data3, aes(x = X, ymin = sem_min, ymax = sem_max), fill = "black", alpha = 0.3) +
  geom_line(data = data3, aes(x = X, y = Average), color = adjustcolor("black", alpha.f = 1), size = 0.5)

# Add individual mouse lines
#for (x in 1:length(MouseIDs2)) {
#  temp <- data.frame(X = data2$X, Average = data2[, MouseIDs2[x]])
#  p <- p + geom_line(data = temp, aes(x = X, y = Average), color = adjustcolor("black", alpha.f = 0.1), size = 1)
#}

# Display the plot and save it as pdf
p <- p + labs(y= y_title, x = x_title) + xlim + ylim + theme_classic() #Print it to save it
ggsave(Plotname, plot = p, width = 4, height = 4)




#------------------------
#Perform the statistics of the Cumulative Occurence using TWO-WAY ANOVA Tukey Post-Hoc

# Combine your dataframes as before
data1 <- data.frame(Value = rep(data1$X, times = ncol(data1[, !names(data1) %in% c("X", "Average", "sem_min", "sem_max")])),
                      Cumulative_Occurrence = unlist(data1[, !names(data1) %in% c("X", "Average", "sem_min", "sem_max")], use.names = FALSE),
                      Group = rep(Group1, times = length(unlist(data1[, !names(data1) %in% c("X", "Average", "sem_min", "sem_max")], use.names = FALSE))))

data2 <- data.frame(Value = rep(data2$X, times = ncol(data2[, !names(data2) %in% c("X", "Average", "sem_min", "sem_max")])),
                     Cumulative_Occurrence = unlist(data2[, !names(data2) %in% c("X", "Average", "sem_min", "sem_max")], use.names = FALSE),
                     Group = rep(Group2, times = length(unlist(data2[, !names(data2) %in% c("X", "Average", "sem_min", "sem_max")], use.names = FALSE))))

data3 <- data.frame(Value = rep(data3$X, times = ncol(data3[, !names(data3) %in% c("X", "Average", "sem_min", "sem_max")])),
                     Cumulative_Occurrence = unlist(data3[, !names(data3) %in% c("X", "Average", "sem_min", "sem_max")], use.names = FALSE),
                        Group = rep(Group3, times = length(unlist(data3[, !names(data3) %in% c("X", "Average", "sem_min", "sem_max")], use.names = FALSE))))

# Combine the three datasets into one
data_combined <- rbind(data1, data2, data3)

# Ensure Group and Velocity are factors, but keep Value numeric
data_combined$Group <- as.factor(data_combined$Group)
data_combined$Value <- as.factor(data_combined$Value)
data_combined$Cumulative_Occurrence <- as.numeric(data_combined$Cumulative_Occurrence)

# Two-way ANOVA with Value as the dependent variable
anova_result_2way <- aov(Cumulative_Occurrence ~ Group * Value, data = data_combined)
tukey_result <- TukeyHSD(anova_result_2way)

#Tukey <- tukey_result$`Group:Cumulative_Occurrence`
#To_keep <- c()
#for (x in 1:nrow(tukey_result$`Group:Velocity`)) {
#  temp <- rownames(Tukey)[x]
#  temp <- as.character(temp)
#  elements <- strsplit(temp, ":")[[1]]
#  if (length(elements) >= 3) {  # Ensure at least two elements exist
#    First_element <- elements[2]
#    First_element <- strsplit(First_element, "-")[[1]][1]
#    Last_element <- elements[length(elements)]
#    if (First_element == Last_element) {
#      To_keep <- c(To_keep, x)
#    }
#  }
#}
#Tukey <- Tukey[To_keep,]
#Tukey <-as.data.frame(Tukey)
#Tukey <- Tukey[Tukey$`p adj` < 0.05,]


write.csv(tukey_result$Group, paste("2WayAnova_TukeyHSD", File, sep = ""))

}

Cumulative_Occurence_Plot(File = "Peak_Speed_Cumulative_Occurrence.csv",
                          Plotname = "Peak Speed Cumulative Occurrence.pdf",
                          x_title = "Peak Speed (mm/s)",
                          y_title = "Cumulative Occurrence",
                          xlim = xlim(0,300),
                          ylim = ylim(0,1),
                          Group1 = "CoffFontaCasp3", 
                          Group1_Col_Index = 2:4,
                          Group2 = "ConFontaCasp3",
                          Group2_Col_Index = 5:12,
                          Group3 = "NotaCasp3",
                          Group3_Col_Index = 13:20)


Cumulative_Occurence_Plot(File = "Duration_Cumulative_Occurrence.csv",
                          Plotname = "Duration Cumulative Occurrence.pdf",
                          x_title = "Bout Duration (s)",
                          y_title = "Cumulative Occurrence",
                          xlim = xlim(0,3),
                          ylim = ylim(0,1),
                          Group1 = "CoffFontaCasp3", 
                          Group1_Col_Index = 2:4,
                          Group2 = "ConFontaCasp3",
                          Group2_Col_Index = 5:12,
                          Group3 = "NotaCasp3",
                          Group3_Col_Index = 13:20)

Cumulative_Occurence_Plot(File = "Peak_Acceleration_Cumulative_Occurrence.csv",
                          Plotname = "Peak Acceleration Cumulative Occurrence.pdf",
                          x_title = "Peak Acceleration (mm/s2)",
                          y_title = "Cumulative Occurrence",
                          xlim = xlim(0,1500),
                          ylim = ylim(0,1),
                          Group1 = "CoffFontaCasp3", 
                          Group1_Col_Index = 2:4,
                          Group2 = "ConFontaCasp3",
                          Group2_Col_Index = 5:12,
                          Group3 = "NotaCasp3",
                          Group3_Col_Index = 13:20)

Cumulative_Occurence_Plot(File = "Peak_Deceleration__Cumulative_Occurrence.csv",
                          Plotname = "Peak Deceleration Cumulative Occurrence.pdf",
                          x_title = "Peak Deceleration (mm/s2)",
                          y_title = "Cumulative Occurrence",
                          xlim = xlim(-1500, 0),
                          ylim = ylim(0,1),
                          Group1 = "CoffFontaCasp3", 
                          Group1_Col_Index = 2:4,
                          Group2 = "ConFontaCasp3",
                          Group2_Col_Index = 5:12,
                          Group3 = "NotaCasp3",
                          Group3_Col_Index = 13:20)

#------------------------------------------------------------------------

