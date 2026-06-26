library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Striatal spectra/Calb1-Cre_DAT-Flpo_taCasp3/CSVfiles")
#Do the analysis for a number of 25 boxes
number_of_box <- 30

#List the csv files within the document
filenames <- list.files(pattern="*.csv")
nfiles <- length(filenames)
nfiles

#Make the dataframe to store the data
df <- data.frame(1:number_of_box)
colnames(df) <- c("Position")



#------------------------------------------------------------------------------------
for (i in 1:nfiles) {
  setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Striatal spectra/Calb1-Cre_DAT-Flpo_taCasp3/CSVfiles")
  Spectra <- as.data.frame(read.csv(filenames[i]))
  
  #Extract the corresponding Noise csv file and denoise all values
  temp <- filenames[i]
  temp <- gsub("_DV", "", temp)
  temp <- gsub("_ML", "", temp)
  temp <- gsub("_MV-DL", "", temp)
  temp <- paste("Noise_", temp, sep = "")
  
  setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Striatal spectra/Calb1-Cre_DAT-Flpo_taCasp3/Noise")
  temp <- as.data.frame(read.csv(temp))

  temp <- Spectra$Value - mean(temp$Mean) #Substract each pixel values by the average noise value
  temp[temp < 0] <- 0 #Replace value below 0 to 0
  Spectra$Value <- temp #Replace the Values by the denoised values
  
  #Crop to get a number of pixels that is mulplicater of 25
  kept_pixels <- floor(length(Spectra[,1])/number_of_box) *number_of_box
  Spectra <- Spectra[1:kept_pixels,]
  
  #Extract and average the values for each boxes
  pixels_per_box <- length(Spectra[,1])/number_of_box
  number_of_positions <- number_of_box + 1
  
  positions <- rep("", times=number_of_positions)  # Create an empty vector to store positions
  
  for (x in 1:number_of_positions) {
    positions[x] <- (x-1) * pixels_per_box
  }
  
  position_A <- positions[1:length(positions)-1]
  
  position_B <- positions[2:length(positions)]
  
  average <- rep("", times=number_of_box)
  
  for (x in 1:number_of_box) {
    average[x] <- mean(
      Spectra[position_A[x]:position_B[x],2]
    )
  }
  
  average <- as.numeric(average) #convert the characters of the previous loop to numeric
  df <- cbind(df, average) #Add the column to the df
  
}


#Remove the Position column
df <- df[,-1]
lol <- t(df)

#Rename the row names
rownames(df) <- paste("Average_box#", as.character(1:number_of_box))

#Add the filenames as the column names
colnames(df) <- filenames

#------------------------------------------------------------------------------------

#Add the metadata to the df
#The file names is a bit messy, so we will need to extract the informations we want
#Extract the channels
Channel <- sapply(strsplit(filenames, "_"), `[`, 1)

#Extract the striatal axis
Striatal_Axis <- sapply(strsplit(filenames, "_"), `[`, 2)

#Extract the Groups
Group <- sapply(strsplit(filenames, "_"), `[`, 3)

#Extract the MouseIDs
MouseID <- sapply(strsplit(filenames, "_"), `[`, 4)

#Extract the Slice#
Slice <- sapply(strsplit(filenames, "_"), `[`, 5)

metadata <- rbind(MouseID, Group, Slice, Channel, Striatal_Axis)
colnames(metadata) <- filenames
metadata <- t(metadata)

df <- t(df)

df <- cbind(metadata, df)

df <- as.data.frame(df)

#------------------------------------------------------------------------------------

df$Group <- gsub('ConFon-taCasp3', 'ConFontaCasp3', df$Group)


#Extract the groups
NotaCasp3 <- df[df$Group == "Control",]
CoffFontaCasp3 <- df[df$Group == "CoffFontaCasp3",]
ConFontaCasp3 <- df[df$Group == "ConFontaCasp3",]

#Extract the bacthes of immuno
Batch1_NotaCasp3 <- rbind(NotaCasp3[NotaCasp3$MouseID == c("F574"),],
                          NotaCasp3[NotaCasp3$MouseID == c("F707"),],
                          NotaCasp3[NotaCasp3$MouseID == c("F795"),])

Batch2_NotaCasp3 <- rbind(NotaCasp3[NotaCasp3$MouseID == c("G0262"),],
                          NotaCasp3[NotaCasp3$MouseID == c("G0215"),],
                          NotaCasp3[NotaCasp3$MouseID == c("G0217"),],
                          NotaCasp3[NotaCasp3$MouseID == c("G0350"),])

Batch1_CoffFontaCasp3 <- rbind(CoffFontaCasp3[CoffFontaCasp3$MouseID == c("F658"),],
                               CoffFontaCasp3[CoffFontaCasp3$MouseID == c("F706"),],
                               CoffFontaCasp3[CoffFontaCasp3$MouseID == c("F708"),]
)

Batch2_CoffFontaCasp3 <- rbind(CoffFontaCasp3[CoffFontaCasp3$MouseID == c("G0218"),],
                               CoffFontaCasp3[CoffFontaCasp3$MouseID == c("G0219"),],
                               CoffFontaCasp3[CoffFontaCasp3$MouseID == c("G0351"),])

#Extract the bacthes of immuno
Batch3_NotaCasp3 <- rbind(NotaCasp3[NotaCasp3$MouseID == c("D696"),],
                          NotaCasp3[NotaCasp3$MouseID == c("E701"),],
                          NotaCasp3[NotaCasp3$MouseID == c("E768"),])

Batch4_NotaCasp3 <- rbind(NotaCasp3[NotaCasp3$MouseID == c("F191"),],
                          NotaCasp3[NotaCasp3$MouseID == c("F272"),],
                          NotaCasp3[NotaCasp3$MouseID == c("F273"),],
                          NotaCasp3[NotaCasp3$MouseID == c("F275"),])

Batch3_ConFontaCasp3 <- rbind(ConFontaCasp3[ConFontaCasp3$MouseID == c("E770"),],
                              ConFontaCasp3[ConFontaCasp3$MouseID == c("E801"),]
)

Batch4_ConFontaCasp3 <- rbind(ConFontaCasp3[ConFontaCasp3$MouseID == c("F273"),],
                              ConFontaCasp3[ConFontaCasp3$MouseID == c("F275"),],
                              ConFontaCasp3[ConFontaCasp3$MouseID == c("F283"),],
                              ConFontaCasp3[ConFontaCasp3$MouseID == c("F284"),]
)

Batch1 <- rbind(Batch1_NotaCasp3, Batch1_CoffFontaCasp3)
Batch2 <- rbind(Batch2_NotaCasp3, Batch2_CoffFontaCasp3)
Batch3 <- rbind(Batch3_NotaCasp3, Batch3_ConFontaCasp3)
Batch4 <- rbind(Batch4_NotaCasp3, Batch4_ConFontaCasp3)

Immuno_Batch <- c(rep("Batch1", times = nrow(Batch1)),
                  rep("Batch2", times = nrow(Batch2)),
                  rep("Batch3", times = nrow(Batch3)),
                  rep("Batch4", times = nrow(Batch4))
)

df <- rbind(Batch1, Batch2, Batch3, Batch4)
df <- cbind(Immuno_Batch, df)

#Replace the df so the fluorescent values are numeric
index <- 6 + number_of_box #index for choosing the rows with the values in the dataframe
temp <- df[,7:index]
temp <- as.matrix(temp)
temp <- as.numeric(temp)
df[,7:index] <- temp

#------------------------------------------------------------------------------------

#Define a function that averages the values from the control groups, so it can be used later for normalization

Average <- c()


average_on_CTL <- function(Batch, Channel, Slice, Striatal_Axis) {
  
  for (x in 1:number_of_box) {
    temp <- df[df$Immuno_Batch == Batch & df$Channel == Channel & df$Group == "Control" & df$Slice == Slice & df$Striatal_Axis == Striatal_Axis,7:index]
    temp <- data.frame(lapply(temp, as.numeric))
    Average[x] <- mean(temp[,x])
  }
  
  Average <- as.numeric(Average)
  
  return(Average)
  
}

#------------------------------------------------------------------------------------

#Get the control averages values for all batches, channels, Slices, and axis
Immuno_Batchs <- unique(df$Immuno_Batch)
Channels <- unique(df$Channel)
Slices <- unique(df$Slice)
Axis <- unique(df$Striatal_Axis)

for (a in Immuno_Batchs) {
  for (b in Channels) {
    for (c in Slices) {
      for (d in Axis) {
        
        #Perform the function
        temp <- average_on_CTL(
          Batch = a,
          Channel = b, 
          Slice = c, 
          Striatal_Axis = d)
        
        #Change the temp vector name by its corresponding control name
        assign(paste("Averaged_Control_Value_", a, "_", b, "_", c, "_", d, sep = ""),temp)
      }
    }
  }
}

#------------------------------------------------------------------------------------
#Generate a normalized_df
normalized_df <- df[,1:6]
temp <- df[,7:index]
colnames(temp) <- paste("Normalized_Box#", as.character(1:number_of_box))
normalized_df <- cbind(normalized_df, temp)

#Replace the values in the normalized_df, dividing the value by the max of its corresponding control average
for (a in Immuno_Batchs) {
  for (b in Channels) {
    for (c in Slices) {
      for (d in Axis) {
        
        temp <- normalized_df[normalized_df$Immuno_Batch == a &
                                normalized_df$Channel == b &
                                normalized_df$Slice == c &
                                normalized_df$Striatal_Axis == d
                              , 
                              7:index]
        
        
        normalized_values <- sweep(temp, STATS = max(get(paste("Averaged_Control_Value_",a, "_",b, "_",c, "_",d, sep = ""))), MARGIN=2, FUN = "/")
        
        
        
        #Replace the values in the df
        normalized_df[normalized_df$Immuno_Batch == a & normalized_df$Channel == b & normalized_df$Slice == c & normalized_df$Striatal_Axis == d, 7:index] <- normalized_values
        
        
      }
    }
  }
}


#------------------------------------------------------------------------------------
#Get the averages curves for the normalized values for all
#Extract the values in a vector
Groups <- unique(normalized_df$Group)

Normalized_Curves_Value <- c()


for (a in Groups) {
  for (b in Channels) {
    for (c in Slices) {
      for (d in Axis) {
        
        #Perform the function
        temp <- normalized_df[normalized_df$Group == a & normalized_df$Channel == b & normalized_df$Slice == c & normalized_df$Striatal_Axis == d,7:index]
        temp <- colMeans(temp)
        temp <- as.data.frame(temp)
        temp <- t(temp)
        temp <- cbind(as.character(a),
                      as.character(c),
                      as.character(b),
                      as.character(d),
                      temp)
        
        Normalized_Curves_Value <- rbind(Normalized_Curves_Value, temp)
        
      }
    }
  }
}

temp <- Normalized_Curves_Value[,5:ncol(Normalized_Curves_Value)]
temp <- as.data.frame(apply(temp, 2, as.numeric))
Normalized_Curves_Value <- Normalized_Curves_Value[,1:4]
Normalized_Curves_Value <- cbind(Normalized_Curves_Value, temp)
colnames(Normalized_Curves_Value) <- colnames(df[,3:ncol(df)])


#------------------------------------------------------------------------------------
#Get the 95% confidence interval 

#Define the function for the IC_95
Confidence_Interval_95 <- function(x) {
  alpha = 0.05
  degrees.freedom = nrow(x) - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  
  IC_95 <- t.score * apply(x, MARGIN = 2, FUN =sd) / sqrt(nrow(x))
}

#Define a function for the SEM
SEM <- function(x) {
  SEM <- apply(x, MARGIN = 2, FUN =sd) / sqrt(nrow(x))
}

#For simplicity of the code, just substitute the IC function by SEM function, so I won't replace everything in the following lines
#Confidence_Interval_95 <- SEM

Normalized_Curves_IC95 <- c()

for (a in Groups) {
  for (b in Channels) {
    for (c in Slices) {
      for (d in Axis) {
        
        #Perform the function
        temp <- normalized_df[normalized_df$Group == a & normalized_df$Channel == b & normalized_df$Slice == c & normalized_df$Striatal_Axis == d,7:index]
        
        temp <- Confidence_Interval_95(temp)
        temp <- as.data.frame(temp)
        temp <- t(temp)
        temp <- cbind(as.character(a),
                      as.character(c),
                      as.character(b),
                      as.character(d),
                      temp)
        
        Normalized_Curves_IC95 <- rbind(Normalized_Curves_IC95, temp)
        
      }
    }
  }
}

temp <- Normalized_Curves_IC95[,5:ncol(Normalized_Curves_IC95)]
temp <- as.data.frame(apply(temp, 2, as.numeric))
Normalized_Curves_IC95 <- Normalized_Curves_IC95[,1:4]
Normalized_Curves_IC95 <- cbind(Normalized_Curves_IC95, temp)
colnames(Normalized_Curves_IC95) <- colnames(df[,3:ncol(df)])

#Generate two dataframes with the maximum and mimimum coordinates of the IC interval
#Starting with ymax
Normalized_Curves_IC95_max <- c()

for (a in Groups) {
  for (b in Channels) {
    for (c in Slices) {
      for (d in Axis) {
        
        #Perform the function
        temp <- Normalized_Curves_Value[Normalized_Curves_Value$Group == a & Normalized_Curves_Value$Channel == b & Normalized_Curves_Value$Slice == c & Normalized_Curves_Value$Striatal_Axis == d, 5:ncol(Normalized_Curves_Value)]
        temp2 <- Normalized_Curves_IC95[Normalized_Curves_IC95$Group == a & Normalized_Curves_IC95$Channel == b & Normalized_Curves_IC95$Slice == c & Normalized_Curves_IC95$Striatal_Axis == d, 5:ncol(Normalized_Curves_IC95)]
        
        temp <- temp + temp2
        temp <- cbind(as.character(a),
                      as.character(c),
                      as.character(b),
                      as.character(d),
                      temp)
        
        Normalized_Curves_IC95_max <- rbind(Normalized_Curves_IC95_max, temp)
        
      }
    }
  }
}

colnames(Normalized_Curves_IC95_max) <- colnames(Normalized_Curves_IC95)

#Then ymin
Normalized_Curves_IC95_min <- c()

for (a in Groups) {
  for (b in Channels) {
    for (c in Slices) {
      for (d in Axis) {
        
        #Perform the function
        temp <- Normalized_Curves_Value[Normalized_Curves_Value$Group == a & Normalized_Curves_Value$Channel == b & Normalized_Curves_Value$Slice == c & Normalized_Curves_Value$Striatal_Axis == d, 5:ncol(Normalized_Curves_Value)]
        temp2 <- Normalized_Curves_IC95[Normalized_Curves_IC95$Group == a & Normalized_Curves_IC95$Channel == b & Normalized_Curves_IC95$Slice == c & Normalized_Curves_IC95$Striatal_Axis == d, 5:ncol(Normalized_Curves_IC95)]
        
        temp <- temp - temp2
        temp <- cbind(as.character(a),
                      as.character(c),
                      as.character(b),
                      as.character(d),
                      temp)
        
        Normalized_Curves_IC95_min <- rbind(Normalized_Curves_IC95_min, temp)
        
      }
    }
  }
}
colnames(Normalized_Curves_IC95_min) <- colnames(Normalized_Curves_IC95)

#Replace the negative values of min values of the IC95 to 0 instead
temp <- Normalized_Curves_IC95_min[,5:ncol(Normalized_Curves_IC95_min)]
temp[temp < 0] <- 0
Normalized_Curves_IC95_min[,5:ncol(Normalized_Curves_IC95_min)] <- temp



#------------------------------------------------------------------------------------

graph <- function(Slice, Channel, Striatal_Axis) {
  


Average <- Normalized_Curves_Value[Normalized_Curves_Value$Group == "Control" &
                              Normalized_Curves_Value$Slice == Slice &
                              Normalized_Curves_Value$Channel == Channel &
                              Normalized_Curves_Value$Striatal_Axis == Striatal_Axis,
                              5:ncol(Normalized_Curves_Value)]

Average2 <- Normalized_Curves_Value[Normalized_Curves_Value$Group == "CoffFontaCasp3" &
                                Normalized_Curves_Value$Slice == Slice &
                                Normalized_Curves_Value$Channel == Channel &
                                Normalized_Curves_Value$Striatal_Axis == Striatal_Axis,
                                5:ncol(Normalized_Curves_Value)]

Average3 <- Normalized_Curves_Value[Normalized_Curves_Value$Group == "ConFontaCasp3" &
                                      Normalized_Curves_Value$Slice == Slice &
                                      Normalized_Curves_Value$Channel == Channel &
                                      Normalized_Curves_Value$Striatal_Axis == Striatal_Axis,
                                    5:ncol(Normalized_Curves_Value)]

Average <- t(Average)
Average <- cbind(Average, c(1:number_of_box), rep("NotaCasp3", times = number_of_box))
colnames(Average) <- c("Normalized_Value", "Position", "Group")

Average2 <- t(Average2)
Average2 <- cbind(Average2, c(1:number_of_box), rep("CoffFon-taCasp3", times = number_of_box))
colnames(Average2) <- c("Normalized_Value", "Position", "Group")

Average3 <- t(Average3)
Average3 <- cbind(Average3, c(1:number_of_box), rep("ConFon-taCasp3", times = number_of_box))
colnames(Average3) <- c("Normalized_Value", "Position", "Group")

Average <- rbind(Average, Average2, Average3)
Average <- as.data.frame(Average)
Average$Normalized_Value <- as.numeric(Average$Normalized_Value)
Average$Position <- as.numeric(Average$Position)




up <- Normalized_Curves_IC95_max[Normalized_Curves_IC95_max$Group == "Control" &
                                    Normalized_Curves_IC95_max$Slice == Slice &
                                    Normalized_Curves_IC95_max$Channel == Channel &
                                    Normalized_Curves_IC95_max$Striatal_Axis == Striatal_Axis,
                               5:ncol(Normalized_Curves_IC95_max)]

up2 <- Normalized_Curves_IC95_max[Normalized_Curves_IC95_max$Group == "CoffFontaCasp3" &
                                   Normalized_Curves_IC95_max$Slice == Slice &
                                   Normalized_Curves_IC95_max$Channel == Channel &
                                   Normalized_Curves_IC95_max$Striatal_Axis == Striatal_Axis,
                                 5:ncol(Normalized_Curves_IC95_max)]

up3 <- Normalized_Curves_IC95_max[Normalized_Curves_IC95_max$Group == "ConFontaCasp3" &
                                    Normalized_Curves_IC95_max$Slice == Slice &
                                    Normalized_Curves_IC95_max$Channel == Channel &
                                    Normalized_Curves_IC95_max$Striatal_Axis == Striatal_Axis,
                                  5:ncol(Normalized_Curves_IC95_max)]


up <- t(up)
up <- cbind(up, c(1:number_of_box), rep("NotaCasp3", times = number_of_box))
colnames(up) <- c("CI_Max", "Position", "Group")

up2 <- t(up2)
up2 <- cbind(up2, c(1:number_of_box), rep("CoffFon-taCasp3", times = number_of_box))
colnames(up2) <- c("CI_Max", "Position", "Group")

up3 <- t(up3)
up3 <- cbind(up3, c(1:number_of_box), rep("ConFon-taCasp3", times = number_of_box))
colnames(up3) <- c("CI_Max", "Position", "Group")

up <- rbind(up, up2, up3)
up <- as.data.frame(up)
up$CI_Max <- as.numeric(up$CI_Max)
up$Position <- as.numeric(up$Position)

low <- Normalized_Curves_IC95_min[Normalized_Curves_IC95_min$Group == "Control" &
                                   Normalized_Curves_IC95_min$Slice == Slice &
                                   Normalized_Curves_IC95_min$Channel == Channel &
                                   Normalized_Curves_IC95_min$Striatal_Axis == Striatal_Axis,
                                 5:ncol(Normalized_Curves_IC95_min)]

low2 <- Normalized_Curves_IC95_min[Normalized_Curves_IC95_min$Group == "CoffFontaCasp3" &
                                    Normalized_Curves_IC95_min$Slice == Slice &
                                    Normalized_Curves_IC95_min$Channel == Channel &
                                    Normalized_Curves_IC95_min$Striatal_Axis == Striatal_Axis,
                                  5:ncol(Normalized_Curves_IC95_min)]

low3 <- Normalized_Curves_IC95_min[Normalized_Curves_IC95_min$Group == "ConFontaCasp3" &
                                     Normalized_Curves_IC95_min$Slice == Slice &
                                     Normalized_Curves_IC95_min$Channel == Channel &
                                     Normalized_Curves_IC95_min$Striatal_Axis == Striatal_Axis,
                                   5:ncol(Normalized_Curves_IC95_min)]

low <- t(low)
low <- cbind(low, c(1:number_of_box), rep("NotaCasp3", times = number_of_box))
colnames(low) <- c("CI_Min", "Position", "Group")

low2 <- t(low2)
low2 <- cbind(low2, c(1:number_of_box), rep("CoffFon-taCasp3", times = number_of_box))
colnames(low2) <- c("CI_Min", "Position", "Group")

low3 <- t(low3)
low3 <- cbind(low3, c(1:number_of_box), rep("ConFon-taCasp3", times = number_of_box))
colnames(low3) <- c("CI_Min", "Position", "Group")

low <- rbind(low, low2, low3)
low <- as.data.frame(low)
low$CI_Min <- as.numeric(low$CI_Min)
low$Position <- as.numeric(low$Position)



data <- data.frame(Average$Group, Average$Position, Average$Normalized_Value, up$CI_Max, low$CI_Min)
colnames(data) <- c("Group", "Position", "Normalized_Value", "CI_Max", "CI_Min")  
  
  
ggplot(data, aes(x = Position, y = Normalized_Value, color = Group)) +
  geom_ribbon(aes(ymin = CI_Min, ymax = CI_Max, fill = Group), alpha = 0.3) +
  xlim(1, number_of_box) + 
  ylim(0,2.5) +
  geom_line() 

}

#------------------------------------------------------------------------------------

graph_individual_mouseID <- function(MouseID, Slice, Channel, Striatal_Axis) {
  
  
  lol <- normalized_df[normalized_df$MouseID == MouseID &
                         normalized_df$Slice == Slice &
                         normalized_df$Channel == Channel &
                         normalized_df$Striatal_Axis == Striatal_Axis,
                       7:ncol(normalized_df)]
  lol <- t(lol)
  lol <- cbind(lol, c(1:number_of_box))
  colnames(lol) <-  c("Normalized_Value", "Position")
  
  ggplot(lol, aes(x = Position, y = Normalized_Value)) +
    xlim(1, number_of_box) + 
    ylim(0,3) 
    geom_line() 
  
  
}

#------------------------------------------------------------------------------------


graph_raw_values <- function(MouseID, Slice, Channel, Striatal_Axis) {
  
  
  lol <- df[df$MouseID == MouseID &
                        df$Slice == Slice &
                        df$Channel == Channel &
                        df$Striatal_Axis == Striatal_Axis,
                       7:ncol(df)]
  lol <- t(lol)
  lol <- cbind(lol, c(1:number_of_box))
  colnames(lol) <-  c("Normalized_Value", "Position")
  
  ggplot(lol, aes(x = Position, y = Normalized_Value)) +
    xlim(1, number_of_box) + 
    ylim(0,1000) +
    geom_line() 
  
  
}

#------------------------------------------------------------------------------------
#Function for figure
screening <- function(MouseID) {
  
  
  a <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice1",
                                Channel = "EYFP",
                                Striatal_Axis = "DV")
  
  b <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice2",
                                Channel = "EYFP",
                                Striatal_Axis = "DV")
  
  c <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice3",
                                Channel = "EYFP",
                                Striatal_Axis = "DV")
  
  d <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice4",
                                Channel = "EYFP",
                                Striatal_Axis = "DV")
  
  e <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice1",
                                Channel = "EYFP",
                                Striatal_Axis = "ML")
  
  f <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice2",
                                Channel = "EYFP",
                                Striatal_Axis = "ML")
  
  g <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice3",
                                Channel = "EYFP",
                                Striatal_Axis = "ML")
  
  h <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice4",
                                Channel = "EYFP",
                                Striatal_Axis = "ML")

  
  i <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice1",
                                Channel = "mCherry",
                                Striatal_Axis = "DV")
  
  j <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice2",
                                Channel = "mCherry",
                                Striatal_Axis = "DV")
  
  k <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice3",
                                Channel = "mCherry",
                                Striatal_Axis = "DV")
  
  l <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice4",
                                Channel = "mCherry",
                                Striatal_Axis = "DV")
  
  m <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice1",
                                Channel = "mCherry",
                                Striatal_Axis = "ML")
  
  n <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice2",
                                Channel = "mCherry",
                                Striatal_Axis = "ML")
  
  o <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice3",
                                Channel = "mCherry",
                                Striatal_Axis = "ML")
  
  p <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice4",
                                Channel = "mCherry",
                                Striatal_Axis = "ML")
  q <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice1",
                                Channel = "TH",
                                Striatal_Axis = "DV")
  
  r <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice2",
                                Channel = "TH",
                                Striatal_Axis = "DV")
  
  s <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice3",
                                Channel = "TH",
                                Striatal_Axis = "DV")
  
  t <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice4",
                                Channel = "TH",
                                Striatal_Axis = "DV")
  
  u <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice1",
                                Channel = "TH",
                                Striatal_Axis = "ML")
  
  v <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice2",
                                Channel = "TH",
                                Striatal_Axis = "ML")
  
  w <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice3",
                                Channel = "TH",
                                Striatal_Axis = "ML")
  
  x <- graph_individual_mouseID(MouseID = MouseID,
                                Slice = "Slice4",
                                Channel = "TH",
                                Striatal_Axis = "ML")  
  
  figure <- ggarrange(a,
                      b,
                      c,
                      d,
                      e,
                      f,
                      g,
                      h,
                      i,
                      j,
                      k,
                      l,
                      m,
                      n,
                      o,
                      p,
                      q,
                      r,
                      s,
                      t,
                      u,
                      v,
                      w,
                      x,
                      ncol = 4, nrow = 6)
  
  
}

#MouseIDs <- unique(df$MouseID)


#for (a in MouseIDs) {
  
  
#png(filename = paste(a, ".png", sep = ""), width = 1000, height = 1000)
  
#plot <- screening(MouseID = a)
#print(plot)
  
#dev.off()
#}
#------------------------------------------------------------------------------------

#Might have to redo the DV axis, but ML axis looks good


Slice_number <- "Slice4"

  
  a <- graph(Slice = Slice_number,
             Channel = "EYFP",
             Striatal_Axis = "DV")
  
  b <- graph(Slice = Slice_number,
             Channel = "mCherry",
             Striatal_Axis = "DV")
  
  c <- graph(Slice = Slice_number,
             Channel = "TH",
             Striatal_Axis = "DV")
  
  
  
  d <- graph(Slice = Slice_number,
             Channel = "EYFP",
             Striatal_Axis = "ML")
  
  e <- graph(Slice = Slice_number,
             Channel = "mCherry",
             Striatal_Axis = "ML")
  
  f <- graph(Slice = Slice_number,
             Channel = "TH",
             Striatal_Axis = "ML")
  

  
  g <- graph(Slice = Slice_number,
             Channel = "EYFP",
             Striatal_Axis = "MV-DL")
  
  h <- graph(Slice = Slice_number,
             Channel = "mCherry",
             Striatal_Axis = "MV-DL")
  
  i <- graph(Slice = Slice_number,
             Channel = "TH",
             Striatal_Axis = "MV-DL")

  

  
  
  
  
  
  
  figure <- ggarrange(a,
                      b,
                      c,
                      d,
                      e,
                      f,
                      g,
                      h,
                      i,
                      ncol = 3, nrow = 3)
figure  
  

#------------------------------------------------------------------------------------
#Perform the statistics of the Cumulative Occurence using TWO-WAY ANOVA Tukey Post-Hoc
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Striatal spectra/Calb1-Cre_DAT-Flpo_taCasp3/Stats")

Slices <- unique(normalized_df$Slice)
Channels <- unique(normalized_df$Channel)
Axis <- unique(normalized_df$Striatal_Axis)

for (a in Slices) {
  for (b in Channels) {
    for (c in Axis) {
      
      
      
      data <- normalized_df[normalized_df$Slice == a & normalized_df$Channel == b & normalized_df$Striatal_Axis == c,]
      
      Bin <- rep(1:number_of_box, times = length(data$Group))
      Group <- rep(data$Group, times = c(rep(number_of_box, times = length(data$Group))))
      Value <- as.vector(as.matrix(t(data[,7:ncol(data)])))
      
      data_combined <- data.frame(Bin, Group, Value)
      
      # Ensure Group and Velocity are factors, but keep Value numeric
      data_combined$Group <- as.factor(data_combined$Group)
      data_combined$Bin <- as.factor(data_combined$Bin)
      data_combined$Value <- as.numeric(data_combined$Value)
      
      # Two-way ANOVA with Value as the dependent variable
      anova_result_2way <- aov(Value ~ Group * Bin, data = data_combined)
      tukey_result <- TukeyHSD(anova_result_2way)
      Tukey <- tukey_result$`Group:Bin`
      To_keep <- c()
      for (x in 1:nrow(Tukey)) {
        temp <- rownames(Tukey)[x]
        temp <- as.character(temp)
        elements <- strsplit(temp, ":")[[1]]
        if (length(elements) >= 3) {  # Ensure at least two elements exist
          First_element <- elements[2]
          First_element <- strsplit(First_element, "-")[[1]][1]
          Last_element <- elements[length(elements)]
          if (First_element == Last_element) {
            To_keep <- c(To_keep, x)
          }
        }
      }
      Tukey <- Tukey[To_keep,]
      Tukey <-as.data.frame(Tukey)
      
      write.csv(Tukey, file = paste("TukeyHSD_perBin_",a, "_", b, "_", c,".csv", sep = ""))
      
      
    }
  }
}
