# Load the packages with library()
library("tidyverse")
library("sf")
library("tmap")
library("nngeo")
library("spdep")
library("sp")
library("data.table")
library('SpatialEpi')
library('INLA')
library("RColorBrewer")
library("grid")

################################ SET DIRECTORY  ########################################
setwd('/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/FINALDFS')

########################## INDEPENDENT VARIABLE DATASET  ###############################
cross_sectional_20222023 <- read.csv("msoaAgg20222023.csv", row.names = NULL, check.names = FALSE)

# drop any row with an NA anywhere
cross_sectional_20222023 <- cross_sectional_20222023 %>% drop_na()
# check length
nrow(cross_sectional_20222023)

# convert rates to decimals
names(cross_sectional_20222023)
cross_sectional_20222023$TotalAchievers <- round(cross_sectional_20222023$TotalAchievers)
cross_sectional_20222023$RateofStudentAchievers <- cross_sectional_20222023$MSOA_GCSE_Mean / 100
cross_sectional_20222023$RateOFFemaleStudents <- cross_sectional_20222023$PercentFemale / 100
cross_sectional_20222023$RateOfStudentReceivingFreeSchoolMeals <- cross_sectional_20222023$PercentOfStudentReceivingFreeSchoolMeals / 100
cross_sectional_20222023$FSMRate <- cross_sectional_20222023$RateOfStudentReceivingFreeSchoolMeals
cross_sectional_20222023$BlackPopulationRate <- cross_sectional_20222023$BlackPopulationRate/100
cross_sectional_20222023$AsianPopulationRate <- cross_sectional_20222023$AsianPopulationRate/100
cross_sectional_20222023$WhitePopulationRate <- cross_sectional_20222023$WhitePopulationRate/100
cross_sectional_20222023$OtherRacePopulationRate <- cross_sectional_20222023$OtherRacePopulationRate/100
cross_sectional_20222023$YPAccessibilityDecile <- cross_sectional_20222023$BordaAccessibilityDecile
cross_sectional_20222023$RateOfStudentUnderAchievers <- (cross_sectional_20222023$TotalPupils_GCSE - cross_sectional_20222023$TotalAchievers) / cross_sectional_20222023$TotalPupils_GCSE

# drop schools that are all female or all male
cross_sectional_20222023 <- cross_sectional_20222023[cross_sectional_20222023$RateOFFemaleStudents > 0.05 & cross_sectional_20222023$RateOFFemaleStudents < 0.95, ]

# drop reference columns not needed for Urbanization and Race
cross_sectional_20222023 <- cross_sectional_20222023[, !names(cross_sectional_20222023) %in% c("LargeCentralMetro", "WhitePopulationRate")]

# remove duplicate MSOA (by MSOA (code))
cross_sectional_20222023 <- cross_sectional_20222023[!duplicated(cross_sectional_20222023$`MSOA (code)`), ]
# check length
nrow(cross_sectional_20222023)

########################## DATA CLEANSING AND WRANGLING  ###############################
# load in the informational data for spatial distribtion maps
spatialDistributionMaps <- read.csv("msoaAgg20222023.csv")
# drop nas again 
spatialDistributionMaps <- spatialDistributionMaps %>% drop_na()
# load in the shapefile
englandMSOA <- st_read("/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/RData/msoa/MSOA_2021_EW_BGC_V3.shp")
englandMSOA <- englandMSOA[!grepl("^W", englandMSOA$MSOA21CD), ]
englandOutline <- st_union(englandMSOA)
plot(englandOutline)

# Data preparation
# estimate the risk of GCSE underperformance from the number of achievers per at risk population 
cross_sectional_20222023$Grade5Underperformance <- cross_sectional_20222023$TotalPupils_GCSE - cross_sectional_20222023$TotalAchievers
cross_sectional_20222023$Expected <- expected(
  population = cross_sectional_20222023$TotalPupils_GCSE,
  cases = cross_sectional_20222023$Grade5Underperformance,
  n.strata = 1
)

# round Expected to a whole number 
cross_sectional_20222023$Expected <- round(cross_sectional_20222023$Expected, .1)

# merge the shapefile and the county level informational data
spatialDistributionMaps <- merge(englandMSOA, spatialDistributionMaps, 
                                 by.x = c("MSOA21CD"),
                                 by.y = c("MSOA..code."))

# merge the shapefile and the cross sectional data
analysis_cross_data <- merge(englandMSOA, cross_sectional_20222023, 
                             by.x = c("MSOA21CD"),
                             by.y = c("MSOA (code)"))

# create an adjacency matrix 
# ensure areas are arranged in order of object id 
analysis_cross_data <- analysis_cross_data[order(analysis_cross_data$MSOA21CD),]
row.names(analysis_cross_data) <- 1:nrow(analysis_cross_data)

# create the neighborhood matrix as an INLA object
adjacencyMatrix <- poly2nb(analysis_cross_data)
# ensure names are carried over into adjacency matrix object
names(adjacencyMatrix) <- analysis_cross_data$MSOA21CD
# check first 10 neighbourhoods, and cross-check with map
head(adjacencyMatrix, n = 10)

# check how many islands - lots of msoas with no schools
sum(sapply(adjacencyMatrix, function(x) identical(x, 0L) || identical(x, 0)))

# create adjacency matrix into a readable graph structure
nb2INLA("adjacencyObject.adj", adjacencyMatrix)
g <- inla.read.graph(filename = "adjacencyObject.adj")

###########################  SPATIAL DISTRIBUTION MAPS ################################## 

# spatial distribution of raw grade 5 achievement rate
sdMSOAAchievementRate <- tm_shape(spatialDistributionMaps) + 
  tm_polygons(
    fill = "MSOA_GCSE_Mean",
    col_alpha = 0.2,
    fill.scale = tm_scale_intervals(n = 5, style = "pretty", values = brewer.pal(5, "RdPu")),
    fill.legend = tm_legend(title = "Rate of Secondary Students Per MSOA Achieving a Grade 5+ on Eng/Math GCSEs", frame = FALSE)
  ) +
  tm_compass(position = c(0.93, 1.2)) +
  tm_scalebar(position = c(0.01, 0.05)) +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.25,
    frame = FALSE)

sdMSOAAchievementRate 

################################# Moran's I ##############################################
# use analysis_cross_data, since that's the dataframe with RateOfStudentUnderAchievers
spatialDistributionMapSP <- as_Spatial(analysis_cross_data, IDs = analysis_cross_data$MSOA21CD)
# inspect
class(spatialDistributionMapSP)
# create an nb object
spatialDistributionMapNB <- poly2nb(spatialDistributionMapSP, row.names = spatialDistributionMapSP$MSOA21CD)
# inspect
class(spatialDistributionMapNB)
# inspect
str(spatialDistributionMapNB, list.len = 10)
# create the list weights object
# allow for islands with no neighborhood weighting
nb_weights_list <- nb2listw(spatialDistributionMapNB, style = 'W', zero.policy = TRUE)
# inspect
class(nb_weights_list)
# Moran's I
mi_value <- moran(spatialDistributionMapSP$RateOfStudentUnderAchievers, nb_weights_list,
                  n = length(nb_weights_list$neighbours), S0 = Szero(nb_weights_list))
# inspect
mi_value
# run a Monte Carlo simulation 599 times
mc_model <- moran.mc(spatialDistributionMapSP$RateOfStudentUnderAchievers, nb_weights_list, nsim = 599)
# inspect
mc_model


##################################### LISA ##############################################
# Local Moran's I
localMoranEngGCSEAchievements <- localmoran(spatialDistributionMapSP$MSOA_GCSE_Mean, nb_weights_list)
# rescale
spatialDistributionMapSP$scaleMSOA_GCSE_Mean <- scale(spatialDistributionMapSP$MSOA_GCSE_Mean)
# create a spatial lag variable 
spatialDistributionMapSP$lag_scale_MSOA_GCSE_Mean <- lag.listw(nb_weights_list, spatialDistributionMapSP$scaleMSOA_GCSE_Mean)
# convert to sf
EngSFmoran_stats <- st_as_sf(spatialDistributionMapSP)
# MAP 2: VERSION WITH SIGNIFICANCE
# set a significance value
sig_level <- 0.1
# classification with significance value
EngSFmoran_stats$quad_sig <- ifelse(EngSFmoran_stats$scaleMSOA_GCSE_Mean > 0 & 
                                      EngSFmoran_stats$lag_scale_MSOA_GCSE_Mean > 0 & 
                                      localMoranEngGCSEAchievements[,5] <= sig_level, 
                                    'high-high', 
                                    ifelse(EngSFmoran_stats$scaleMSOA_GCSE_Mean <= 0 & 
                                             EngSFmoran_stats$lag_scale_MSOA_GCSE_Mean <= 0 & 
                                             localMoranEngGCSEAchievements[,5] <= sig_level, 
                                           'low-low', 
                                           ifelse(EngSFmoran_stats$scaleMSOA_GCSE_Mean > 0 & 
                                                    EngSFmoran_stats$lag_scale_MSOA_GCSE_Mean <= 0 & 
                                                    localMoranEngGCSEAchievements[,5] <= sig_level, 
                                                  'high-low', 
                                                  ifelse(EngSFmoran_stats$scaleMSOA_GCSE_Mean <= 0 & 
                                                           EngSFmoran_stats$lag_scale_MSOA_GCSE_Mean > 0 & 
                                                           localMoranEngGCSEAchievements[,5] <= sig_level, 
                                                         'low-high',
                                                         'not-significant'))))

# map only the statistically significant results here
LISAMap <- tm_shape(EngSFmoran_stats) +
  tm_polygons(fill = "quad_sig", col_alpha = 0.2,
              fill.scale = tm_scale_categorical(values = c("#de2d26", "#fee0d2", "#deebf7","#3182bd","white")),
              fill.legend = tm_legend(frame = FALSE, title = "Clusters")) +
  tm_compass(position = c(0.93, 1.25)) +
  tm_scalebar(position = c(0.01, 0.05)) +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.25,
    frame = FALSE)

LISAMap

################################# CORRELATION MATRIX #######################################
cor_matrix <- cor(cross_sectional_20222023[, c("BordaAccessibilityDecile", "FSMRate", "IDACI", "SAMHI", "RateOFFemaleStudents",
                                               "LargeFringeMetro", "MediumMetro", "SmallMetro", "Micropolitan", "Noncore", "BlackPopulationRate", "AsianPopulationRate", "OtherRacePopulationRate")],
                  use = "complete.obs")
round(cor_matrix, 2)
library(corrplot)
corrplot(cor_matrix, method = "color", type = "full", tl.cex = 0.7)

# drop IDACI due to multicolinearity 
cross_sectional_20222023 <- cross_sectional_20222023[, !names(cross_sectional_20222023) %in% c("IDACI")]
analysis_cross_data  <- analysis_cross_data[, !names(analysis_cross_data) %in% c("IDACI")]

################################# BAYESIAN MODEL SET UP #######################################
# unstructured and structured data
analysis_cross_data$id_area_structured <- 1:nrow(analysis_cross_data)
analysis_cross_data$id_area_unstructured <- 1:nrow(analysis_cross_data)

# regresssion model 
unperformanceRisk <- Grade5Underperformance ~ 1 + YPAccessibilityDecile + FSMRate + SAMHI + RateOFFemaleStudents +
  LargeFringeMetro + MediumMetro + SmallMetro + Micropolitan + Noncore + BlackPopulationRate + AsianPopulationRate + OtherRacePopulationRate +
  f(id_area_structured, model = "besag", graph = g, scale.model = TRUE) + f(id_area_unstructured, model = "iid")

# fit the model and store the results
results_cross_sectional <- inla(unperformanceRisk, 
                                family = "poisson", 
                                data = analysis_cross_data, 
                                E = Expected,
                                control.predictor = list(compute = TRUE), 
                                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE))

# print the results 
summary(results_cross_sectional)

################################# REPORTING RISKS #######################################
# calculate the regression coefficients as relative risk
# extract the results to report them as relative risks
options(scipen = 7)
fixed <- results_cross_sectional$summary.fixed
RR_table <- exp(fixed[, c("mean", "0.025quant", "0.5quant", "0.975quant")])
# use rows 2, 3 and 4
RR_table[c(2:14),]

################## POSTERIOR DISTRIBUTIONS AND EXCEEDANCE PROBABILITY ####################
names(results_cross_sectional$marginals.fixed)
beta_YPAccessibilityDecile <- inla.smarginal(results_cross_sectional$marginals.fixed$YPAccessibilityDecile)
beta_FSMRate <- inla.smarginal(results_cross_sectional$marginals.fixed$FSMRate)
beta_SAMHI <- inla.smarginal(results_cross_sectional$marginals.fixed$SAMHI)
beta_RateOFFemaleStudents <- inla.smarginal(results_cross_sectional$marginals.fixed$RateOFFemaleStudents)
beta_LargeFringeMetro <- inla.smarginal(results_cross_sectional$marginals.fixed$LargeFringeMetro)
beta_MediumMetro <- inla.smarginal(results_cross_sectional$marginals.fixed$MediumMetro)
beta_SmallMetro <- inla.smarginal(results_cross_sectional$marginals.fixed$SmallMetro)
beta_Micropolitan <- inla.smarginal(results_cross_sectional$marginals.fixed$Micropolitan)
beta_Noncore <- inla.smarginal(results_cross_sectional$marginals.fixed$Noncore)
beta_BlackPopulationRate <- inla.smarginal(results_cross_sectional$marginals.fixed$BlackPopulationRate)
beta_AsianPopulationRate <- inla.smarginal(results_cross_sectional$marginals.fixed$AsianPopulationRate)
beta_OtherRacePopulationRate <- inla.smarginal(results_cross_sectional$marginals.fixed$OtherRacePopulationRate)

post_dist_beta <- data.frame(beta_YPAccessibilityDecile, beta_FSMRate, beta_SAMHI, beta_RateOFFemaleStudents,
                             beta_LargeFringeMetro, beta_MediumMetro, beta_SmallMetro,
                             beta_Micropolitan, beta_Noncore, beta_BlackPopulationRate, beta_AsianPopulationRate,
                             beta_OtherRacePopulationRate)

colnames(post_dist_beta)[1]  <- "beta_YPAccessibilityDecile_coef"
colnames(post_dist_beta)[2]  <- "beta_YPAccessibilityDecile_density"
colnames(post_dist_beta)[3]  <- "beta_FSMRate_coef"
colnames(post_dist_beta)[4]  <- "beta_FSMRate_density"
colnames(post_dist_beta)[5]  <- "beta_SAMHI_coef"
colnames(post_dist_beta)[6]  <- "beta_SAMHI_density"
colnames(post_dist_beta)[7]  <- "beta_RateOFFemaleStudents_coef"
colnames(post_dist_beta)[8]  <- "beta_RateOFFemaleStudents_density"
colnames(post_dist_beta)[9]  <- "beta_LargeFringeMetro_coef"
colnames(post_dist_beta)[10] <- "beta_LargeFringeMetro_density"
colnames(post_dist_beta)[11] <- "beta_MediumMetro_coef"
colnames(post_dist_beta)[12] <- "beta_MediumMetro_density"
colnames(post_dist_beta)[13] <- "beta_SmallMetro_coef"
colnames(post_dist_beta)[14] <- "beta_SmallMetro_density"
colnames(post_dist_beta)[15] <- "beta_Micropolitan_coef"
colnames(post_dist_beta)[16] <- "beta_Micropolitan_density"
colnames(post_dist_beta)[17] <- "beta_Noncore_coef"
colnames(post_dist_beta)[18] <- "beta_Noncore_density"
colnames(post_dist_beta)[19] <- "beta_BlackPopulationRate_coef"
colnames(post_dist_beta)[20] <- "beta_BlackPopulationRate_density"
colnames(post_dist_beta)[21] <- "beta_AsianPopulationRate_coef"
colnames(post_dist_beta)[22] <- "beta_AsianPopulationRate_density"
colnames(post_dist_beta)[23] <- "beta_OtherRacePopulationRate_coef"
colnames(post_dist_beta)[24] <- "beta_OtherRacePopulationRate_density"

post_dist_beta$beta_YPAccessibilityDecile_coef <- exp(post_dist_beta$beta_YPAccessibilityDecile_coef)
post_dist_beta$beta_FSMRate_coef <- exp(post_dist_beta$beta_FSMRate_coef)
post_dist_beta$beta_SAMHI_coef <- exp(post_dist_beta$beta_SAMHI_coef)
post_dist_beta$beta_RateOFFemaleStudents_coef <- exp(post_dist_beta$beta_RateOFFemaleStudents_coef)
post_dist_beta$beta_LargeFringeMetro_coef <- exp(post_dist_beta$beta_LargeFringeMetro_coef)
post_dist_beta$beta_MediumMetro_coef <- exp(post_dist_beta$beta_MediumMetro_coef)
post_dist_beta$beta_SmallMetro_coef <- exp(post_dist_beta$beta_SmallMetro_coef)
post_dist_beta$beta_Micropolitan_coef <- exp(post_dist_beta$beta_Micropolitan_coef)
post_dist_beta$beta_Noncore_coef <- exp(post_dist_beta$beta_Noncore_coef)
post_dist_beta$beta_BlackPopulationRate_coef <- exp(post_dist_beta$beta_BlackPopulationRate_coef)
post_dist_beta$beta_AsianPopulationRate_coef <- exp(post_dist_beta$beta_AsianPopulationRate_coef)
post_dist_beta$beta_OtherRacePopulationRate_coef <- exp(post_dist_beta$beta_OtherRacePopulationRate_coef)

################# YP Accessibility Decile ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$YPAccessibilityDecile)
x <- post_dist_beta$beta_YPAccessibilityDecile_coef
y <- post_dist_beta$beta_YPAccessibilityDecile_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 1.0093221); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 1.0093221, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 1.0059692, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0126864, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_YP <- 1 - inla.pmarginal(1, marg); prob_YP

################# FSM Rate ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$FSMRate)
x <- post_dist_beta$beta_FSMRate_coef
y <- post_dist_beta$beta_FSMRate_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 2.6983735); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 2.6983735, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 2.5072541, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 2.9041349, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_FSM <- 1 - inla.pmarginal(1, marg); prob_FSM

################# Percent Female ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$RateOFFemaleStudents)
x <- post_dist_beta$beta_RateOFFemaleStudents_coef
y <- post_dist_beta$beta_RateOFFemaleStudents_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 0.7210111); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 0.7210111, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.6251018, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 0.8315707, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_Female <- 1 - inla.pmarginal(1, marg); prob_Female

################# SAMHI ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$SAMHI)
x <- post_dist_beta$beta_SAMHI_coef
y <- post_dist_beta$beta_SAMHI_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 1.0130215); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 1.0130215, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 1.0086818, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0173818, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_SAMHI <- 1 - inla.pmarginal(1, marg); prob_SAMHI

################# Large Fringe Metro ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$LargeFringeMetro)
x <- post_dist_beta$beta_LargeFringeMetro_coef
y <- post_dist_beta$beta_LargeFringeMetro_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 1.0629991); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 1.0629991, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 1.0294274, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0976689, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_LFM <- 1 - inla.pmarginal(1, marg); prob_LFM

################# Medium Metro ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$MediumMetro)
x <- post_dist_beta$beta_MediumMetro_coef
y <- post_dist_beta$beta_MediumMetro_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 0.9399271); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 0.9399271, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.9064616, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 0.9746279, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_MM <- 1 - inla.pmarginal(1, marg); prob_MM

################# Small Metro ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$SmallMetro)
x <- post_dist_beta$beta_SmallMetro_coef
y <- post_dist_beta$beta_SmallMetro_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 1.0297664); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 1.0297664, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.9828944, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0788780, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_SM <- 1 - inla.pmarginal(1, marg); prob_SM

################# Micropolitan ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$Micropolitan)
x <- post_dist_beta$beta_Micropolitan_coef
y <- post_dist_beta$beta_Micropolitan_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 0.9872682); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 0.9872682, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.9313690, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0465590, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_MI <- 1 - inla.pmarginal(1, marg); prob_MI

################# Noncore ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$Noncore)
x <- post_dist_beta$beta_Noncore_coef
y <- post_dist_beta$beta_Noncore_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 1.0171285); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 1.0171285, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.9453317, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0943886, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_NC <- 1 - inla.pmarginal(1, marg); prob_NC

################# BlackPopulationRate ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$BlackPopulationRate)
x <- post_dist_beta$beta_BlackPopulationRate_coef
y <- post_dist_beta$beta_BlackPopulationRate_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 1.3362929); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 1.3362929, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 1.0707915, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.6675256, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_BP <- 1 - inla.pmarginal(1, marg); prob_BP

################# AsianPopulationRate ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$AsianPopulationRate)
x <- post_dist_beta$beta_AsianPopulationRate_coef
y <- post_dist_beta$beta_AsianPopulationRate_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 0.9955803); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 0.9955803, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.9226101, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.0743378, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_AP <- 1 - inla.pmarginal(1, marg); prob_AP

################# OtherRacePopulationRate ####################
marg <- inla.tmarginal(exp, results_cross_sectional$marginals.fixed$OtherRacePopulationRate)
x <- post_dist_beta$beta_OtherRacePopulationRate_coef
y <- post_dist_beta$beta_OtherRacePopulationRate_density
plot(x, y, type = "l", xlab = "Estimated Relative Risk", ylab = "Posterior Probability Density (Plausibility)")
idx <- which(x >= 0.1489270); polygon(c(x[idx], rev(x[idx])), c(y[idx], rep(0, length(idx))), col = "#fee0d2", border = NA)
abline(v = 0.1489270, lty = "dashed", col = "darkgrey", lwd = 4)
abline(v = 0.1013889, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 0.2187450, lty = "dashed", col = "darkgrey", lwd = 2)
abline(v = 1.00, col = "red", lwd = 2)
prob_OP <- 1 - inla.pmarginal(1, marg); prob_OP

############################# PROBABILITIES ###############################
prob_YP
prob_FSM
prob_Female
prob_SAMHI
prob_LFM
prob_MM
prob_SM
prob_MI
prob_NC
prob_BP
prob_AP
prob_OP
RR_table[c(2:14),]

############################# MAPPING THE RELATIVE RISKS ###############################

# extract area specific risks and store them in riskratio
riskratio <- results_cross_sectional$summary.fitted.values
head(riskratio, n=10)

# slot the area-specific risk estimates into the spatial data frame object analysis_cross_data
analysis_cross_data$RR <- riskratio[, "mean"]       # Relative risk
analysis_cross_data$LL <- riskratio[, "0.025quant"] # Lower credibility limit
analysis_cross_data$UL <- riskratio[, "0.975quant"] # Upper credibility limit

# check
head(analysis_cross_data)

# distribution for area-specific relative risks to inform how we create the risk categories for our map and legends
summary(analysis_cross_data$RR)

# create buckets 
RiskCategorylist <- c("<0.25", "0.26 to 0.5", "0.51 to 0.75", "0.76 to 0.99", "1.00-1.009 (null value)",
                      "1.01 to 1.25", "1.26 to 1.50", "1.51 to 1.75", "1.76 to 2.0", "2.01+")

# color palette scheme 
RRPalette <- c("darkblue","#33a6fe","#98cffe","#cbe6fe","#fef9f9","#fed5d5","#feb1b1","#fe8e8e","#fe2424", "darkred")

analysis_cross_data$RelativeRiskCat <- NA
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 0 & analysis_cross_data$RR <= 0.25] <- -4
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 0.25 & analysis_cross_data$RR <= 0.50] <- -3
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 0.50 & analysis_cross_data$RR <= 0.75] <- -2
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 0.75 & analysis_cross_data$RR < 1] <- -1
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR>= 1 & analysis_cross_data$RR < 1.01] <- 0
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR>= 1.01 & analysis_cross_data$RR <= 1.25] <- 1
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 1.25 & analysis_cross_data$RR <= 1.50] <- 2
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 1.50 & analysis_cross_data$RR <= 1.75] <- 3
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 1.75 & analysis_cross_data$RR <= 2] <- 4
analysis_cross_data$RelativeRiskCat[analysis_cross_data$RR> 2 & analysis_cross_data$RR <= 3] <- 5


analysis_cross_data$Significance <- NA
analysis_cross_data$Significance[analysis_cross_data$LL<1 & analysis_cross_data$UL>1] <- 0    # NOT SIGNIFICANT
analysis_cross_data$Significance[analysis_cross_data$LL==1 | analysis_cross_data$UL==1] <- 0  # NOT SIGNIFICANT
analysis_cross_data$Significance[analysis_cross_data$LL>1 & analysis_cross_data$UL>1] <- 1    # SIGNIFICANT INCREASE
analysis_cross_data$Significance[analysis_cross_data$LL<1 & analysis_cross_data$UL<1] <- -1   # SIGNIFICANT DECREASE

# computes exceedance probabilities
analysis_cross_data$ExceedProb <- sapply(results_cross_sectional$marginals.fitted.values,
                                         FUN = function(marg){1 - inla.pmarginal(q = 1.00, marginal = marg)})
# rounds it to three decimal places
analysis_cross_data$ExceedProb <- round(analysis_cross_data$ExceedProb, 3)

# categorisation
analysis_cross_data$ProbCat <- NA
analysis_cross_data$ProbCat[analysis_cross_data$ExceedProb>=0 & analysis_cross_data$ExceedProb< 0.01] <- 1
analysis_cross_data$ProbCat[analysis_cross_data$ExceedProb>=0.01 & analysis_cross_data$ExceedProb< 0.20] <- 2
analysis_cross_data$ProbCat[analysis_cross_data$ExceedProb>=0.20 & analysis_cross_data$ExceedProb< 0.40] <- 3
analysis_cross_data$ProbCat[analysis_cross_data$ExceedProb>=0.40 & analysis_cross_data$ExceedProb< 0.60] <- 4
analysis_cross_data$ProbCat[analysis_cross_data$ExceedProb>=0.60 & analysis_cross_data$ExceedProb< 0.80] <- 5
analysis_cross_data$ProbCat[analysis_cross_data$ExceedProb>=0.80 & analysis_cross_data$ExceedProb<= 1.00] <- 6

# labeling for legend
ProbCategorylist <- c("<0.01", "0.01-0.20", "0.20-0.39", "0.40-0.59", "0.60-0.79", "0.80-1.00")

gm <- c("Bolton","Bury","Manchester","Oldham","Rochdale",
        "Salford","Stockport","Tameside","Trafford","Wigan")
pattern <- paste0("^(", paste(gm, collapse = "|"), ") [0-9]")
manchester <- englandMSOA[grepl(pattern, englandMSOA$MSOA21NM), ]

nrow(manchester)

boroughs <- c("City of London","Barking and Dagenham","Barnet","Bexley","Brent","Bromley","Camden","Croydon","Ealing","Enfield","Greenwich","Hackney","Hammersmith and Fulham","Haringey","Harrow","Havering","Hillingdon","Hounslow","Islington","Kensington and Chelsea","Kingston upon Thames","Lambeth","Lewisham","Merton","Newham","Redbridge","Richmond upon Thames","Southwark","Sutton","Tower Hamlets","Waltham Forest","Wandsworth","Westminster")
pattern <- paste0("^(", paste(boroughs, collapse = "|"), ") [0-9]")
london <- englandMSOA[grepl(pattern, englandMSOA$MSOA21NM), ]

nrow(london)

manchester_rr <- analysis_cross_data[analysis_cross_data$MSOA21CD %in% manchester$MSOA21CD, ]
london_rr <- analysis_cross_data[analysis_cross_data$MSOA21CD %in% london$MSOA21CD, ]

nrow(manchester_rr)
nrow(london_rr)

############################# Map A ############################# 
# England 
sdRelativeRisk <- tm_shape(englandOutline) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.01) +
  tm_shape(analysis_cross_data) +
  tm_polygons(
    fill = "RelativeRiskCat",
    col  = "grey60",
    lwd  = .1,
    fill.scale = tm_scale_categorical(
      values = RRPalette,
      labels = RiskCategorylist
    ),
    fill.legend = tm_legend(title = "Relative Risk: GCSE Underperformance", frame = FALSE)
  ) +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(
    legend.outside = TRUE, legend.outside.position = "right",
    legend.outside.size = 0.25, legend.title.size = 0.7,
    legend.text.size = 0.6, frame = FALSE
  )
sdRelativeRisk
tmap_save(sdRelativeRisk, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/RelativeRisk/england20222023.png", width = 2400, height = 2400, dpi = 300)

# manchester
map_manchester <- tm_shape(manchester) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.1) +
  tm_shape(manchester_rr) +
  tm_polygons(
    fill = "RelativeRiskCat",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(values = RRPalette, labels = RiskCategorylist),
    fill.legend = tm_legend_hide()
  ) +
  tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
  tm_layout(frame = TRUE, bg.color = "white")

# london 
map_london <- tm_shape(london) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.1) +
  tm_shape(london_rr) +
  tm_polygons(
    fill = "RelativeRiskCat",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(values = RRPalette, labels = RiskCategorylist),
    fill.legend = tm_legend_hide()
  ) +
  tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
  tm_layout(frame = TRUE, bg.color = "white")

map_manchester
tmap_save(map_manchester, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/RelativeRisk/manchester20222023.png", width = 2400, height = 2400, dpi = 300)

map_london
tmap_save(map_london, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/RelativeRisk/london20222023.png", width = 2400, height = 2400, dpi = 300)

###################################################################
map_A <- tm_shape(analysis_cross_data) + 
  tm_polygons(fill = "RelativeRiskCat", 
              fill.scale = tm_scale_categorical(values = RRPalette, labels = RiskCategorylist),
              fill.legend = tm_legend(frame = FALSE, "Relative Risk: GCSE Underperformance")) +
  tm_compass(position = c(0.93, 1.15)) +
  tm_scalebar(position = c(0.01, 0.05)) +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.25,
    frame = FALSE)

map_A

############################## MAP B #####################################
map_B <- tm_shape(englandOutline) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.01) +
  tm_shape(analysis_cross_data) +
  tm_polygons(
    fill = "Significance",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(
      values = c("blue", "white", "red"),
      labels = c("Decreased Risk: Significant", "Not Significant", "Increased Risk: Significant")
    ),
    fill.legend = tm_legend(title = "Significance", frame = FALSE)
  ) +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(
    legend.outside = TRUE, legend.outside.position = "right",
    legend.outside.size = 0.25, legend.title.size = 0.7,
    legend.text.size = 0.6, frame = FALSE
  )
map_B
tmap_save(map_B, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/Significance/england20222023.png", width = 2400, height = 2400, dpi = 300)

manchester_sig <- analysis_cross_data[analysis_cross_data$MSOA21CD %in% manchester$MSOA21CD, ]
london_sig <- analysis_cross_data[analysis_cross_data$MSOA21CD %in% london$MSOA21CD, ]

map_B_manchester <- tm_shape(manchester) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.1) +
  tm_shape(manchester_sig) +
  tm_polygons(
    fill = "Significance",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(
      values = c("blue", "white", "red"),
      labels = c("Decreased Risk: Significant", "Not Significant", "Increased Risk: Significant")
    ),
    fill.legend = tm_legend_hide()
  ) +
  tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
  tm_layout(frame = TRUE, bg.color = "white")

map_B_london <- tm_shape(london) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.1) +
  tm_shape(london_sig) +
  tm_polygons(
    fill = "Significance",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(
      values = c("blue", "white", "red"),
      labels = c("Decreased Risk: Significant", "Not Significant", "Increased Risk: Significant")
    ),
    fill.legend = tm_legend_hide()
  ) +
  tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
  tm_layout(frame = TRUE, bg.color = "white")

map_B_manchester
tmap_save(map_B_manchester, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/Significance/manchester20222023.png", width = 2400, height = 2400, dpi = 300)
map_B_london
tmap_save(map_B_london, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/Significance/london20222023.png", width = 2400, height = 2400, dpi = 300)


map_B <- tm_shape(analysis_cross_data) + 
  tm_polygons(fill = "Significance", 
              fill.scale = tm_scale_categorical(values = c("blue", "white", "red"), 
                                                labels = c("Decreased Risk: Significant", "Not Significant", "Increased Risk: Significant")),
              fill.legend = tm_legend(frame = FALSE, title = "Significance")) +
  tm_compass(position = c(0.93, 1.15)) +
  tm_scalebar(position = c(0.01, 0.05)) +
  tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.25,
    frame = FALSE)

map_B

############################## MAP C #####################################

map_C <- tm_shape(englandOutline) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.01) +
  tm_shape(analysis_cross_data) +
  tm_polygons(
    fill = "ProbCat",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(values = "brewer.reds", labels = ProbCategorylist),
    fill.legend = tm_legend(title = "Exceedance Probability: P(RR > 1.00)", frame = FALSE)
  ) +
  tm_compass(type = "arrow", position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(
    legend.outside = TRUE, legend.outside.position = "right",
    legend.outside.size = 0.25, legend.title.size = 0.7,
    legend.text.size = 0.6, frame = FALSE
  )
map_C
tmap_save(map_C, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/ExceedanceProbability/england20222023.png", width = 2400, height = 2400, dpi = 300)

manchester_prob <- analysis_cross_data[analysis_cross_data$MSOA21CD %in% manchester$MSOA21CD, ]
london_prob <- analysis_cross_data[analysis_cross_data$MSOA21CD %in% london$MSOA21CD, ]

map_C_manchester <- tm_shape(manchester) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.1) +
  tm_shape(manchester_prob) +
  tm_polygons(
    fill = "ProbCat",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(values = "brewer.reds", labels = ProbCategorylist),
    fill.legend = tm_legend_hide()
  ) +
  tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
  tm_layout(frame = TRUE, bg.color = "white")

map_C_london <- tm_shape(london) +
  tm_polygons(fill = "grey85", col = "grey60", lwd = 0.1) +
  tm_shape(london_prob) +
  tm_polygons(
    fill = "ProbCat",
    col  = "grey60", lwd = .1,
    fill.scale = tm_scale_categorical(values = "brewer.reds", labels = ProbCategorylist),
    fill.legend = tm_legend_hide()
  ) +
  tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
  tm_layout(frame = TRUE, bg.color = "white")

map_C_manchester
tmap_save(map_C_manchester, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/ExceedanceProbability/manchester20222023.png", width = 2400, height = 2400, dpi = 300)
map_C_london
tmap_save(map_C_london, filename = "/Users/maggiebickerstaffe/Desktop/Dissertation/WalkableCharities/BayesianRisks/ExceedanceProbability/london20222023.png", width = 2400, height = 2400, dpi = 300)

tmap_arrange(map_A, map_B, map_C, nrow = 1)

############################# EXPORT RESULTS TO CSV ###############################

# numeric codes
analysis_cross_data$RiskCategoryLabelMapA <- factor(
  analysis_cross_data$RelativeRiskCat,
  levels = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
  labels = RiskCategorylist
)

analysis_cross_data$SignificanceLabelMapB <- factor(
  analysis_cross_data$Significance,
  levels = c(-1, 0, 1),
  labels = c("Significant Decrease", "Not Significant", "Significant Increase")
)

analysis_cross_data$ProbCategoryLabelMapC <- factor(
  analysis_cross_data$ProbCat,
  levels = c(1, 2, 3, 4, 5, 6),
  labels = ProbCategorylist
)

# drop geometry to get a plain data frame for export
export_df <- st_drop_geometry(analysis_cross_data)

# write to csv
write.csv(export_df, "relative_risk_results_20222023.csv", row.names = FALSE)

################################# SAVE RDS #######################################
saveRDS(results_cross_sectional, "results_cross_sectional_2022_2023.rds")

################################# VIF #######################################
library(car)
vif_model <- lm(Grade5Underperformance ~ YPAccessibilityDecile + FSMRate + PercentFemale +
                  SAMHI + LargeFringeMetro + MediumMetro + SmallMetro + Micropolitan + Noncore,
                data = cross_sectional_20222023)
vif(vif_model)

################################# MODEL COMPARISON #######################################
hist(results_cross_sectional$cpo$pit,
     breaks = 20,
     main = "PIT Histogram (Model Calibration Check)",
     xlab = "PIT value",
     col = "lightblue",
     border = "white")

fixed_pred <- results_cross_sectional$summary.linear.predictor[, "mean"]

fixed_only <- inla(Grade5Underperformance ~ YPAccessibilityDecile + FSMRate + PercentFemale +
                     SAMHI + LargeFringeMetro + MediumMetro + SmallMetro + Micropolitan + Noncore,
                   family = "poisson", E = cross_sectional_20222023$Expected, data = cross_sectional_20222023,
                   control.compute = list(dic = TRUE, waic = TRUE))

fitted_fixed <- fixed_only$summary.fitted.values[, "mean"] * cross_sectional_20222023$Expected
cor(cross_sectional_20222023$Grade5Underperformance, fitted_fixed)^2

