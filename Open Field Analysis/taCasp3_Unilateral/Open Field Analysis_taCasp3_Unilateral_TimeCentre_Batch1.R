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

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1")
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
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1/TimeCentre_Output")

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
Perc_Centre$Group <- sapply(strsplit(names(DLC_Traces), "_Bilateral_"), "[", 1)
Perc_Centre$MouseID <- sapply(strsplit(names(DLC_Traces), "_Bilateral_"), "[", 2)
rownames(Perc_Centre) <- names(DLC_Traces)

for (MouseID in names(DLC_Traces)) {
  Centre <- DLC_Traces[[MouseID]] %>%
    dplyr::filter(body_x > left_corner & body_x < right_corner) %>%
    dplyr::filter(body_y > bottom_corner & body_y < top_corner)
  
  Perc_Centre[MouseID,]$Perc_Centre <- (nrow(Centre)/nrow(DLC_Traces[[MouseID]]))*100
}

#Export to csv
write.csv(Perc_Centre, "taCasp3Bilateral_Batch1_OpenField_%Centre.csv")

#----------------------------------------------------------------------------------
#Plot
setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/OpenField_Velocity Analysis/taCasp3_Batch/Unilateral/Batch1/TimeCentre_Output")
Toplot <- rbind(read.csv("taCasp3Bilateral_Batch1_OpenField_%Centre.csv"),
                read.csv("taCasp3Bilateral_Batch2_OpenField_%Centre.csv"),
                read.csv("taCasp3Bilateral_Batch3_OpenField_%Centre.csv"))
Toplot$MouseID <- sapply(strsplit(Toplot$Group, "_Unilateral_"), "[", 2)
Toplot$MouseID <- sapply(strsplit(Toplot$MouseID, ".csv"), "[", 1)
Toplot$Group <- sapply(strsplit(Toplot$Group, "_Unilateral_"), "[", 1)

Toplot <- Toplot %>%
  mutate(Sex = case_when(
    MouseID == "F574" ~ "M",
    MouseID == "F707" ~ "M",
    MouseID == "F658" ~ "M",
    MouseID == "F708" ~ "M",
    MouseID == "F574" ~ "M",
    MouseID == "F706" ~ "F",
    MouseID == "G0262" ~ "M",
    MouseID == "G0215" ~ "F",
    MouseID == "G0217" ~ "F",
    MouseID == "G0350" ~ "F",
    MouseID == "G0264" ~ "M",
    MouseID == "G0218" ~ "F",
    MouseID == "G0219" ~ "F",
    MouseID == "G0351" ~ "F",
    MouseID == "I0135" ~ "M",
    MouseID == "I0137" ~ "M",
    MouseID == "I0214" ~ "M",
    MouseID == "I0029" ~ "F",
    MouseID == "I0031" ~ "F",
    MouseID == "I0250" ~ "F",
    MouseID == "I0138" ~ "M",
    MouseID == "I0212" ~ "M",
    MouseID == "I0213" ~ "M",
    MouseID == "I0032" ~ "F",
    MouseID == "I0248" ~ "F",
    MouseID == "I0256" ~ "F"
  ))

Toplot$Group <- factor(Toplot$Group, levels = c("NotaCasp3", "CoffFontaCasp3"))
Toplot |>
  tidyplot(x = Group, y = Perc_Centre, color = Group) |>
  adjust_colors(c("NotaCasp3" = "#B2B2B2", "CoffFontaCasp3" = "#FDB3B3")) |> 
  add_boxplot() |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  adjust_y_axis_title("% Time in Centre") |>
  add_test_asterisks(hide_info = TRUE, method = "t.test") |> #Perform T-test
  adjust_x_axis_title("")
