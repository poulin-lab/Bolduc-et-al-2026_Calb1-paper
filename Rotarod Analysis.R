library(ggplot2)
library(tidyplots)
library(lme4)
library(emmeans)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(WRS2)

#library(lmerTest) # Pour obtenir les p-values

setwd("C:/Users/cyril/OneDrive - McGill University (1)/Documents/PhD/Raw data/Rotarod_Slope")

#-------------------------------------------------------------------------------------------
############################################################################################
############################################################################################
###################################     taCasp3    #########################################
############################################################################################
############################################################################################

df <- read.csv(file = "taCasp3(1).csv")
#Generate empty dataframe to store the data
Slope_df <- data.frame(matrix(nrow = length(unique(df$MouseID)) * length(unique(df$Day)),
                              ncol = 5))

colnames(Slope_df) <- c("Day", "Group", "MouseID", "Sex", "Slope")
Slope_df$Day <- rep(unique(df$Day), times = length(unique(df$MouseID)))
Slope_df$MouseID <- rep(unique(df$MouseID), each = length(unique(df$Day)))

#Extract the intraday slope
for (Mouse_n in unique(df$MouseID)) {
  for (Day_n in unique(df$Day)) {
    Sample <- df %>%
      filter(Day == Day_n, Mouse_n == MouseID)
    
    Trials <- Sample$Trial
    Performance <- Sample$Latency
    model <- lm(Performance ~ Trials)
    Slope_df[Slope_df$Day == Day_n & Slope_df$MouseID == Mouse_n,]$Slope <- coef(model)[2]
    
    
  }
}

#Add a column of the groups and sex to the Slope_df
mouse_map <- df |> 
  dplyr::select(MouseID, Group, Sex) |> 
  unique()

Slope_df <- Slope_df |>
  dplyr::select(-Group, -Sex) |> 
  left_join(mouse_map, by = "MouseID")


#Plot
Slope_df$Group <- factor(Slope_df$Group, levels = c("NotaCasp3", "CreOFFFlpoON-taCasp3", "CreONFlpoON-taCasp3"))
df$Group <- factor(df$Group, levels = c("NotaCasp3", "CreOFFFlpoON-taCasp3", "CreONFlpoON-taCasp3"))

#Plot initial coordination
Toplot <- df %>% dplyr::filter(Trial %in% c(1,6,11))
Toplot$Trial <- as.character(Toplot$Trial)
Toplot <- Toplot %>%mutate(Trial = recode(Trial, 
                        "1" = "Day1",
                        "6" = "Day2",
                        "11" = "Day3"))
Toplot |>
  tidyplot(x = Trial, y = Latency, color = Group) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1) |>
  adjust_y_axis_title("Initial coordination (sec.)") |>
  adjust_x_axis_title("") |>
  add_boxplot() 

Toplot |>
  tidyplot(x = Group, y = Latency, color = Trial) |>
  # Dans tidyplot, on passe le vecteur de couleurs directement ici :
  adjust_colors(c("Day1" = "#B2B2B2", "Day2" = "#FDB3B3", "Day3" = "#B2DCB4")) |> 
  add_boxplot() |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  adjust_y_axis_title("Initial coordination (sec.)") |>
  adjust_x_axis_title("")

#Plot training slope per day 
Toplot <- Slope_df
Toplot$Day <- as.character(Toplot$Day)
Toplot |>
  tidyplot(x = Day, y = Slope, color = Group) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1) |>
  adjust_y_axis_title("Learning slope (sec./trial)") |>
  add_boxplot() 


Toplot |>
  tidyplot(x = Group, y = Slope, color = Day) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16, color = "black") |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1, color = "black") |>
  adjust_y_axis_title("Learning slope (sec./trial)") |>
  add_boxplot() 

#-------------------------------------------------------------------------------------------
#Perform ANOVA
############################################################################################
##################################  FOR THE INI. CORD.  ####################################
############################################################################################
df <- df %>%
  mutate(Group = recode(Group, 
                        "NotaCasp3" = "NotaCasp3",
                        "CreOFFFlpoON-taCasp3" = "CoffFontaCasp3",
                        "CreONFlpoON-taCasp3" = "ConFontaCasp3"))

IniCord <- df %>% dplyr::filter(Trial %in% c(1,6,11))
IniCord$Trial <- as.character(IniCord$Trial)
IniCord <- IniCord %>%mutate(Trial = recode(Trial, 
                                            "1" = "Day1",
                                            "6" = "Day2",
                                            "11" = "Day3"))

#Check assumption
#Outliers
IniCord %>%
  group_by(Trial,Group) %>%
  identify_outliers(Latency)
#Note: No extreme outliers were found

#Normality
ggqqplot(IniCord, "Latency", ggtheme = theme_bw()) +
  facet_grid(Group ~ Trial)
#Note: All the points fall approximately along the reference line, for each cell. So we can assume normality of the data.

#Homogeneity of variance assumption
IniCord %>%
  group_by(Trial)%>%
  levene_test(Latency ~ Group)
#The Levene’s test is significant, therefore assuming there the homogeneity criteria is unmet,
#so I will perform a robust 2way anova

#Robust statiscal test
IniCord$Group <- as.factor(IniCord$Group)
IniCord$Trial <- as.factor(IniCord$Trial)

#By trial
res.robust <- bwtrim(Latency ~ Group * Trial, id = MouseID, data = IniCord, tr = 0.2)
res.robust
#Note: There was a significant day and group effect and interaction effect
#Since the homogeneity criteria was unmet, I will use a GamesHowell test

GH_Test <- IniCord %>%
  group_by(Trial) %>%
  games_howell_test(Latency ~ Group)

GH_Test

#Note: The mixed two ANOVA shows no significant pvalues
write.csv(GH_Test, "taCasp3_InitialCoordination_WelchAOV_GamesHowell.csv")


#By group
res.robust <- bwtrim(Latency ~ Trial * Group, id = MouseID, data = IniCord, tr = 0.2)
res.robust
#Note: There was a significant day and group effect and interaction effect
#Since the homogeneity criteria was unmet, I will use a GamesHowell test

GH_Test <- IniCord %>%
  group_by(Group) %>%
  games_howell_test(Latency ~ Trial)

GH_Test

#Note: The mixed two ANOVA shows no significant pvalues
write.csv(GH_Test, "taCasp3_InitialCoordination_WelchAOV_GamesHowell_ByGroup.csv")

############################################################################################
##################################  FOR THE SLOPE  #########################################
############################################################################################

Slope_df <- Slope_df %>%
  mutate(Group = recode(Group, 
                        "NotaCasp3" = "NotaCasp3",
                        "CreOFFFlpoON-taCasp3" = "CoffFontaCasp3",
                        "CreONFlpoON-taCasp3" = "ConFontaCasp3"))

#Check assumption
#Outliers
Slope_df %>%
  group_by(Day,Group) %>%
  identify_outliers(Slope)


Slope_df <- Slope_df %>% dplyr::filter(!MouseID == "M12") #Remove "M12" since it is an extreme outlier

#Normality
Slope_df %>%
  group_by(Day,Group) %>%
  shapiro_test(Slope)

ggqqplot(Slope_df, "Slope", ggtheme = theme_bw()) +
  facet_grid(Group ~ Day)
#Note: The data are normally distributed saphiro wilk indicates a pvalue of near 0 for Day CoffFon-taCasp3
#that seems to be explained by the fact that the data are really close from each other by visual inspection.
#Based on that visual inspection of QQplot, I consider the data normally distributed

#Homogeneity of variance assumption
Slope_df %>%
  group_by(Day)%>%
  levene_test(Slope ~ Group)
#The Levene’s test is not significant (p > 0.05). Therefore, we can assume the homogeneity of variances in the different groups.

#Homogeneity of covariances assumption
#Compute Box’s M-test:
box_m(Slope_df[, "Slope", drop = FALSE], Slope_df$Group)
#There was homogeneity of covariances, as assessed by Box’s test of equality of covariance matrices (p > 0.001).

# Compute ANOVA
#By day
res.aov <- Slope_df %>% 
  anova_test(dv = Slope, wid = MouseID, within = Day, between = Group)
res.aov
#Note: The mixed two ANOVA shows significant pvalues, so will proceed to Tukey pairwise multiple comparison

pwc_tukey <- Slope_df %>%
  group_by(Day) %>%
  tukey_hsd(Slope ~ Group)

write.csv(pwc_tukey, file = "taCasp3LearningSlope_PerDay_MixedAOV_TukeyHSD.csv")

#By group
Slope_df$Group <- as.factor(Slope_df$Group)
Slope_df$Day <- as.factor(Slope_df$Day)

res.aov <- Slope_df %>% 
  anova_test(dv = Slope, wid = MouseID, within = Day, between = Group)
res.aov
#Note: The mixed two ANOVA shows significant pvalues, so will proceed to Tukey pairwise multiple comparison

pwc_tukey <- Slope_df %>%
  group_by(Group) %>%
  tukey_hsd(Slope ~ Day)
pwc_tukey

write.csv(pwc_tukey, file = "taCasp3LearningSlope_PerDay_MixedAOV_TukeyHSD_byGroup.csv")

#-------------------------------------------------------------------------------------------
############################################################################################
############################################################################################
#################################     taCasp3 Uni   ########################################
############################################################################################
############################################################################################

df <- read.csv(file = "taCasp3_Uni.csv")
#Generate empty dataframe to store the data
Slope_df <- data.frame(matrix(nrow = length(unique(df$MouseID)) * length(unique(df$Day)),
                              ncol = 5))

colnames(Slope_df) <- c("Day", "Group", "MouseID", "Sex", "Slope")
Slope_df$Day <- rep(unique(df$Day), times = length(unique(df$MouseID)))
Slope_df$MouseID <- rep(unique(df$MouseID), each = length(unique(df$Day)))

#Extract the intraday slope
for (Mouse_n in unique(df$MouseID)) {
  for (Day_n in unique(df$Day)) {
    Sample <- df %>%
      filter(Day == Day_n, Mouse_n == MouseID)
    
    Trials <- Sample$Trial
    Performance <- Sample$Latency
    model <- lm(Performance ~ Trials)
    Slope_df[Slope_df$Day == Day_n & Slope_df$MouseID == Mouse_n,]$Slope <- coef(model)[2]
    
    
  }
}

#Add a column of the groups and sex to the Slope_df
mouse_map <- df |> 
  dplyr::select(MouseID, Group, Sex) |> 
  unique()

Slope_df <- Slope_df |>
  dplyr::select(-Group, -Sex) |> 
  left_join(mouse_map, by = "MouseID")


#Plot
Slope_df$Group <- factor(Slope_df$Group, levels = c("NotaCasp3", "CoffFontaCasp3"))
df$Group <- factor(df$Group, levels = c("NotaCasp3", "CoffFontaCasp3"))

#Plot initial coordination
Toplot <- df %>% dplyr::filter(Trial %in% c(1,6,11))
Toplot <- IniCord
Toplot$Trial <- as.character(Toplot$Trial)
Toplot <- Toplot %>%mutate(Trial = recode(Trial, 
                                          "1" = "Day1",
                                          "6" = "Day2",
                                          "11" = "Day3"))
Toplot |>
  tidyplot(x = Trial, y = Latency, color = Group) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1) |>
  adjust_y_axis_title("Initial coordination (sec.)") |>
  adjust_x_axis_title("") |>
  add_boxplot() 

#Plot training slope per day 
Toplot <- Slope_df
Toplot$Day <- as.character(Toplot$Day)
Toplot |>
  tidyplot(x = Day, y = Slope, color = Group) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1) |>
  adjust_y_axis_title("Learning slope (sec./trial)") |>
  add_boxplot() 

#-------------------------------------------------------------------------------------------
#Perform ANOVA
############################################################################################
##################################  FOR THE INI. CORD.  ####################################
############################################################################################

IniCord <- df %>% dplyr::filter(Trial %in% c(1,6,11))
IniCord$Trial <- as.character(IniCord$Trial)
IniCord <- IniCord %>%mutate(Trial = recode(Trial, 
                                            "1" = "Day1",
                                            "6" = "Day2",
                                            "11" = "Day3"))

#Check assumption
#Outliers
IniCord %>%
  group_by(Trial,Group) %>%
  identify_outliers(Latency)
IniCord <- IniCord %>% dplyr::filter(!MouseID %in% c("M2", "M11")) #Remove "M2" and "M11" since they are extreme outlier

#Normality
ggqqplot(IniCord, "Latency", ggtheme = theme_bw()) +
  facet_grid(Group ~ Trial)
#Note: All the points fall approximately along the reference line, for each cell. So we can assume normality of the data.

#Homogeneity of variance assumption
IniCord %>%
  group_by(Trial)%>%
  levene_test(Latency ~ Group)
#The Levene’s test is not significant, so can proceed to 2way anova

# Compute ANOVA
#Bytrial
res.aov <- IniCord %>% 
  anova_test(dv = Latency, wid = MouseID, within = Day, between = Group)
res.aov
#Note: The mixed two ANOVA shows significant pvalues, so will proceed to Tukey pairwise multiple comparison

Tukey <- IniCord %>%
  group_by(Day) %>%
  tukey_hsd(Latency ~ Group)

write.csv(Tukey, "taCasp3Uni_InitialCoordination_2wayMixAOV_Tukey.csv")

#Bygroup
res.aov <- IniCord %>% 
  anova_test(dv = Latency, wid = MouseID, within = Day, between = Group)
res.aov
#Note: The mixed two ANOVA shows significant pvalues, so will proceed to Tukey pairwise multiple comparison
IniCord$Day <- as.factor(IniCord$Day)
IniCord$Group <- as.factor(IniCord$Group)


Tukey <- IniCord %>%
  group_by(Group) %>%
  tukey_hsd(Latency ~ Day)
Tukey

write.csv(Tukey, "taCasp3Uni_InitialCoordination_2wayMixAOV_Tukey_ByGroup.csv")



############################################################################################
##################################  FOR THE SLOPE  #########################################
############################################################################################

#Check assumption
#Outliers
Slope_df %>%
  group_by(Day,Group) %>%
  identify_outliers(Slope)

Slope_df <- Slope_df %>% dplyr::filter(!MouseID == "M11") #Remove "M11" since it is an extreme outlier

#Normality
Slope_df %>%
  group_by(Day,Group) %>%
  shapiro_test(Slope)

ggqqplot(Slope_df, "Slope", ggtheme = theme_bw()) +
  facet_grid(Group ~ Day)
#Note: The data seem normally distributed

#Homogeneity of variance assumption
Slope_df %>%
  group_by(Day)%>%
  levene_test(Slope ~ Group)
#The Levene’s test is not significant (p > 0.05). Therefore, we can assume the homogeneity of variances in the different groups.

# Compute ANOVA
res.aov <- Slope_df %>% 
  anova_test(dv = Slope, wid = MouseID, within = Day, between = Group)
res.aov
#Note: The mixed two ANOVA shows significant pvalues, so will proceed to Tukey pairwise multiple comparison
#By day
Slope_df$Day <- as.factor(Slope_df$Day)
Slope_df$Group <- as.factor(Slope_df$Group)

pwc_tukey <- Slope_df %>%
  group_by(Day) %>%
  tukey_hsd(Slope ~ Group)

write.csv(pwc_tukey, file = "taCasp3_Uni_LearningSlope_PerDay_MixedAOV_TukeyHSD.csv")

#By group
pwc_tukey <- Slope_df %>%
  group_by(Group) %>%
  tukey_hsd(Slope ~ Day)
pwc_tukey

write.csv(pwc_tukey, file = "taCasp3_Uni_LearningSlope_PerDay_MixedAOV_TukeyHSD_ByGroup.csv")

#-------------------------------------------------------------------------------------------
############################################################################################
############################################################################################
###################################     hM4Di    ###########################################
############################################################################################
############################################################################################

df <- read.csv(file = "DREADD(1).csv")
#Generate empty dataframe to store the data
Slope_df <- data.frame(matrix(nrow = length(unique(df$MouseID)) * length(unique(df$Day)),
                              ncol = 5))

colnames(Slope_df) <- c("Day", "Group", "MouseID", "Sex", "Slope")
Slope_df$Day <- rep(unique(df$Day), times = length(unique(df$MouseID)))
Slope_df$MouseID <- rep(unique(df$MouseID), each = length(unique(df$Day)))

#Extract the intraday slope
for (Mouse_n in unique(df$MouseID)) {
  for (Day_n in unique(df$Day)) {
    Sample <- df %>%
      filter(Day == Day_n, Mouse_n == MouseID)
    
    Trials <- Sample$Trial
    Performance <- Sample$Latency
    model <- lm(Performance ~ Trials)
    Slope_df[Slope_df$Day == Day_n & Slope_df$MouseID == Mouse_n,]$Slope <- coef(model)[2]
    
    
  }
}

#Add a column of the groups and sex to the Slope_df
mouse_map <- df |> 
  dplyr::select(MouseID, Group, Sex) |> 
  unique()

Slope_df <- Slope_df |>
  dplyr::select(-Group, -Sex) |> 
  left_join(mouse_map, by = "MouseID")

#Plot
Slope_df$Group <- factor(Slope_df$Group, levels = c("fDIO-mCherry", "CoffFonhM4Di", "ConFonhM4Di"))
df$Group <- factor(df$Group, levels = c("fDIO-mCherry", "CoffFonhM4Di", "ConFonhM4Di"))

#Plot initial coordination
Toplot <- df %>% dplyr::filter(Trial %in% c(1,6,11))
Toplot$Trial <- as.character(Toplot$Trial)
Toplot <- Toplot %>%mutate(Trial = recode(Trial, 
                                          "1" = "Day1",
                                          "6" = "Day2",
                                          "11" = "Day3"))
Toplot |>
  tidyplot(x = Trial, y = Latency, color = Group) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1) |>
  adjust_y_axis_title("Initial coordination (sec.)") |>
  add_boxplot() 

#Plot training slope per day 
Toplot <- Slope_df
Toplot$Day <- as.character(Toplot$Day)
Toplot |>
  tidyplot(x = Day, y = Slope, color = Group) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "M"), shape = 16) |>
  add_data_points_beeswarm(data = filter_rows(Sex == "F"), shape = 1) |>
  adjust_y_axis(limits = c(-27,70)) |>
  adjust_y_axis_title("Learning slope (s / trial)") |>
  add_boxplot() 

#-------------------------------------------------------------------------------------------
#Perform ANOVA
############################################################################################
##################################  FOR THE INI. CORD.  ####################################
############################################################################################

IniCord <- df %>% dplyr::filter(Trial %in% c(1,6,11))
IniCord$Trial <- as.character(IniCord$Trial)
IniCord <- IniCord %>%mutate(Trial = recode(Trial, 
                                          "1" = "Day1",
                                          "6" = "Day2",
                                          "11" = "Day3"))

#Check assumption
#Outliers
IniCord %>%
  group_by(Trial,Group) %>%
  identify_outliers(Latency)


IniCord <- IniCord %>% dplyr::filter(!MouseID %in% c("M10")) #Remove "M10" since it is an extreme outlier

#Normality
ggqqplot(IniCord, "Latency", ggtheme = theme_bw()) +
  facet_grid(Group ~ Trial)
#Note: All the points fall approximately along the reference line, for each cell. So we can assume normality of the data.

#Homogeneity of variance assumption
IniCord %>%
  group_by(Trial)%>%
  levene_test(Latency ~ Group)
#The Levene’s test is significant, therefore assuming there the homogeneity criteria is unmet,
#so I will perform a robust 2way anova

#Robust statiscal test
IniCord$Group <- as.factor(IniCord$Group)
IniCord$Trial <- as.factor(IniCord$Trial)

res.robust <- bwtrim(Latency ~ Group * Trial, id = MouseID, data = IniCord, tr = 0.2)
res.robust
#Note: There was a significant day and group effect but subsignificant interaction effect
#Since the homogeneity criteria was unmet, I will use a GamesHowell test

#By trial
GH_Test <- IniCord %>%
  group_by(Trial) %>%
  games_howell_test(Latency ~ Group)

GH_Test

#Note: The mixed two ANOVA shows no significant pvalues
write.csv(GH_Test, "DREADD_InitialCoordination_WelchAOV_GamesHowell.csv")


#By group
GH_Test <- IniCord %>%
  group_by(Group) %>%
  games_howell_test(Latency ~ Trial)

GH_Test

#Note: The mixed two ANOVA shows no significant pvalues
write.csv(GH_Test, "DREADD_InitialCoordination_WelchAOV_GamesHowell_Bygroup.csv")
############################################################################################
##################################  FOR THE SLOPE  #########################################
############################################################################################

Slope_df <- Slope_df %>%
  mutate(Group = recode(Group, 
                        "NotaCasp3" = "NotaCasp3",
                        "CreOFFFlpoON-taCasp3" = "CoffFontaCasp3",
                        "CreONFlpoON-taCasp3" = "ConFontaCasp3"))

#Check assumption
#Outliers
Slope_df %>%
  group_by(Day,Group) %>%
  identify_outliers(Slope)


Slope_df <- Slope_df %>% dplyr::filter(!MouseID %in% c("M21","M23")) #Remove "M21" and "M23" since it is an extreme outlier

#Normality
Slope_df %>%
  group_by(Day,Group) %>%
  shapiro_test(Slope)

ggqqplot(Slope_df, "Slope", ggtheme = theme_bw()) +
  facet_grid(Group ~ Day)
#Note: All data are normally distributed, except day1 for CoffFon-taCasp3, but the ggqqplot() indicates that it is within the range

#Homogeneity of variance assumption
Slope_df %>%
  group_by(Day)%>%
  levene_test(Slope ~ Group)
#The Levene’s test is not significant (p > 0.05). Therefore, we can assume the homogeneity of variances in the different groups.

#Homogeneity of covariances assumption
#Compute Box’s M-test:
box_m(Slope_df[, "Slope", drop = FALSE], Slope_df$Group)
#There was homogeneity of covariances, as assessed by Box’s test of equality of covariance matrices (p > 0.001).

# Compute ANOVA
res.aov <- Slope_df %>% 
  anova_test(dv = Slope, wid = MouseID, within = Day, between = Group)
res.aov
#Note: The mixed two ANOVA shows no significant pvalues



