library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(DescTools)
library(magick)
library(imager)
library(lattice)
library(viridis)
library(imager)
#library(paletteer)
library(gplots)
library(sommer)
library(png)
library(grid)
library(tidyplots)
library(ggpubr)
library(rstatix)

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/DREADD/Batch1")
#------------------------------------------------------------------------------

#Load the image
image_filenames <- list.files(pattern="*.png")
img <- readPNG(image_filenames[1])
h <- nrow(img)
w <- ncol(img)
windows(width = 10, height = 10 * (h/w)) 
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, type = "n", xlim = c(0, w), ylim = c(h, 0), asp = 1, axes = FALSE)
rasterImage(img, 0, h, w, 0)
#Place three dots on top left and top right corners and bottom left corners respectively
pts <- locator(3)

#Show the points
points(pts$x, pts$y, col = "red", pch = 19, cex = 1.5)

#-----------------------------------------------------------------------------

#Read the CSV files
#List the csv files within the document
filenames <- list.files(pattern="*.csv")
nfiles <- length(filenames)

filenames

#-----------------------------------------------------------------------------
#Import the data
DLC_Traces <- list()

for (MouseID in filenames) {
  #Read the csv
  DLC_Traces[[MouseID]] <- read.csv(file = MouseID)
  
  #Reorder the dataframe
  rownames(DLC_Traces[[MouseID]]) <- DLC_Traces[[MouseID]]$scorer
  colnames(DLC_Traces[[MouseID]]) <- paste(DLC_Traces[[MouseID]]["bodyparts",], DLC_Traces[[MouseID]]["coords",], sep = "_")
  DLC_Traces[[MouseID]] <- DLC_Traces[[MouseID]][-1:-2,]
  DLC_Traces[[MouseID]] <- DLC_Traces[[MouseID]][,-1]
  
  #Filter for high confidence frames
  DLC_Traces[[MouseID]] <- DLC_Traces[[MouseID]] %>%
    dplyr::filter(body_likelihood > 0.9)
  
  #Convert data to numeric
  DLC_Traces[[MouseID]][] <- lapply(DLC_Traces[[MouseID]], function(x) as.numeric(as.character(x)))
  
}

#-----------------------------------------------------------------------------
#Extract the time in centre
deltax <- pts$x[2] - pts$x[1]
deltay <- pts$y[3] - pts$y[1]

#Define the threshold of the corners
n_squares <- 7 #number of squares to define the sides and center of open field
left_corner <- pts$x[1] + deltax/n_squares
right_corner <- pts$x[2] - deltax/n_squares
top_corner <- pts$y[3] - deltay/n_squares
bottom_corner <- pts$y[1] + deltay/n_squares

#Export the time in centre per mouse
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/DREADD/Batch1/TimeCentre_Output")

for (MouseID in names(DLC_Traces)) {
  p <- ggplot(data = DLC_Traces[[MouseID]], aes(x = body_x, y = body_y)) +
    geom_path() +
    annotate(
      "rect",
      xmin = pts$x[1],
      xmax = pts$x[2],
      ymin = pts$y[3],
      ymax = pts$y[1],
      fill = "transparent",
      color = "blue", 
      size = 2) +
    annotate(
      "rect",
      xmin = left_corner,
      xmax = right_corner,
      ymin = bottom_corner,
      ymax = top_corner,
      fill = "transparent",
      color = "blue", 
      size = 2,
      linetype = "dotted") +
    theme_void()
  
  ggsave(paste(MouseID, "Timeincentre.png", sep = "_"), 
         plot = p, 
         bg = "transparent", 
         width = 5, # Optional: specify dimensions
         height = 5)
  
}


#-----------------------------------------------------------------------------
#Extract the time in centre  
Perc_Centre <- data.frame(matrix(nrow = length(DLC_Traces), ncol = 4))
colnames(Perc_Centre) <- c("MouseID", "Sex", "Group", "Perc_Centre")
Perc_Centre$Group <- sapply(strsplit(names(DLC_Traces), "_"), "[", 1)
Perc_Centre$MouseID <- sapply(strsplit(names(DLC_Traces), "_"), "[", 2)
rownames(Perc_Centre) <- names(DLC_Traces)

for (MouseID in names(DLC_Traces)) {
  Centre <- DLC_Traces[[MouseID]] %>%
    dplyr::filter(body_x > left_corner & body_x < right_corner) %>%
    dplyr::filter(body_y > bottom_corner & body_y < top_corner)
  
  Perc_Centre[MouseID,]$Perc_Centre <- (nrow(Centre)/nrow(DLC_Traces[[MouseID]]))*100
}

#Specify the sex
Perc_Centre <- Perc_Centre %>%
  mutate(Sex = case_when(
    MouseID == "I0251.csv" ~ "M",
    MouseID == "I0258.csv" ~ "M",
    MouseID == "I0540.csv" ~ "M",
    MouseID == "I0603.csv" ~ "M",
    MouseID == "I0617.csv" ~ "F",
    MouseID == "I0413.csv" ~ "F",
    MouseID == "I0252.csv" ~ "M",
    MouseID == "I0260.csv" ~ "M",
    MouseID == "I0544.csv" ~ "M",
    MouseID == "I0604.csv" ~ "M",
    MouseID == "I0418.csv" ~ "M",
    MouseID == "I0804.csv" ~ "F",
    MouseID == "I0253.csv" ~ "M",
    MouseID == "I0263.csv" ~ "M",
    MouseID == "I0606.csv" ~ "M",
    MouseID == "I0258.csv" ~ "M",
    MouseID == "I0540.csv" ~ "M",
    MouseID == "I0603.csv" ~ "M",
    MouseID == "I0617.csv" ~ "F",
    MouseID == "I0413.csv" ~ "F",
    MouseID == "I0408.csv" ~ "F",
    MouseID == "I0619.csv" ~ "F"
    
  ))

Perc_Centre

#Export to csv
write.csv(Perc_Centre, "DREADD_Batch1_OpenField_%Centre.csv")

#----------------------------------------------------------------------------------
#Plot
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/DREADD/Batch1/TimeCentre_Output")
Toplot <- rbind(read.csv("DREADD_Batch1_OpenField_%Centre.csv"),
                read.csv("DREADD_Batch2_OpenField_%Centre.csv"))
unique(Toplot$Group)
Toplot$Group <- factor(Toplot$Group, levels = c("fDIO-mCherry", "CoffFon-DREADD", "ConFon-DREADD"))

Toplot |>
  tidyplot(x = Group, y = Perc_Centre, color = Group) |>
  adjust_colors(c("fDIO-mCherry" = "#B2B2B2", "CoffFon-DREADD" = "#FDB3B3", "ConFon-DREADD" = "#B2DCB4")) |> 
  add_boxplot() |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  adjust_y_axis_title("% Time in Centre") |>
  adjust_x_axis_title("")

#----------------------------------------------------------------------------------
#Perform one-way anova

#Summary statictic
Toplot %>%
  group_by(Group) %>%
  get_summary_stats(Perc_Centre, type = "mean_sd")

#Test for outliers
Toplot %>%
  group_by(Group) %>%
  identify_outliers(Perc_Centre)
Toplot <-  Toplot %>%
  dplyr::filter(!MouseID == "J0013.csv") #J0013 is an extreme outlier so will be removed

#Test for normality
model  <- lm(Perc_Centre ~ Group, data = Toplot) # Build the linear model
ggqqplot(residuals(model)) # Create a QQ plot of residuals
shapiro_test(residuals(model)) # Compute saphiro-wilk
#Normality can be assumed

#Test for normality per group
Toplot %>%
  group_by(Group) %>%
  shapiro_test(Perc_Centre)

ggqqplot(Toplot, "Perc_Centre", facet.by = "Group")
#The pvalue and qqplot intragroup is not significant so normality can further be assumed

#Test for homogeneity of variance
Toplot %>% levene_test(Perc_Centre ~ Group)
#pValue < 0.05, so homogeneity of variance criteria is unmet. So, I will compute a Welch one-way ANOVA

#Welch one-way ANOVA
welchAOV <- Toplot %>% welch_anova_test(Perc_Centre ~ Group)
welchAOV

#Note: The pValue was significant, wo will proceed to a Games-Howell post-hoc test
welch_test <- Toplot %>% games_howell_test(Perc_Centre ~ Group)

write.csv(welch_test, "WelchAOV_GamesHowellposthoc_DREADD.csv")
