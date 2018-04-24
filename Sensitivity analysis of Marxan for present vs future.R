##############################################################################################################################
##############################################################################################################################
###Analysing Marxan results and finding top 10%

#
##
###
####
#####
#Clear workspace
rm(list=ls())

require(rgdal)
require(rgeos)
require(tidyverse)
setwd("K:/Christopher_PhD/Github/ParticleTracking") #Need this to access subcode properly

###################Initialize run with these important parameters
my_resolution <- 10000 #defines raster cell size and controls for biased larvae release
#source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")

Depth_class <- "Nearshore"
Baseline_pld <- 56

Top_percent <- 0.5

Depth_class_map <- readOGR(paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/", Depth_class), Depth_class)

########################################################################
########################################################################
###[1] Get overlaps

#Current
Marxan_present <- read.table(paste0("./Marxan/From_R/output/hexagon/", Depth_class, "/pld", Baseline_pld, "/Mar_out_ssoln.dat"))
Marxan_present <- dplyr::rename(Marxan_present, Poly_ID = V1, SSOLN = V2)
Marxan_present <- sp::merge(Depth_class_map, Marxan_present, by = "Poly_ID", all.x = TRUE)
Marxan_present <- Marxan_present@data
Marxan_present$pld <- Baseline_pld
Marxan_present$included <- 1
Marxan_present$RAND <- sample(1:nrow(Marxan_present), nrow(Marxan_present), replace = FALSE)
Marxan_present <- Marxan_present[with(Marxan_present, order(-SSOLN, RAND)),]
Marxan_present <- Marxan_present[c(1:round(Top_percent*nrow(Marxan_present))),]

#Future
Marxan_future <- read.table(paste0("./Marxan_future/From_R/output/hexagon/", Depth_class, "/pld", Baseline_pld, "/Mar_out_ssoln.dat"))
Marxan_future <- dplyr::rename(Marxan_future, Poly_ID = V1, SSOLN = V2)
Marxan_future <- sp::merge(Depth_class_map, Marxan_future, by = "Poly_ID", all.x = TRUE)
Marxan_future <- Marxan_future@data
Marxan_future$pld <- Baseline_pld
Marxan_future$included <- 1
Marxan_future$RAND <- sample(1:nrow(Marxan_future), nrow(Marxan_future), replace = FALSE)
Marxan_future <- Marxan_future[with(Marxan_future, order(-SSOLN, RAND)),]
Marxan_future <- Marxan_future[c(1:round(Top_percent*nrow(Marxan_future))),]

########################################################################
########################################################################
###Top 10% crossover between different plds

present_future_compare <- merge(Marxan_present, Marxan_future, by = "Poly_ID")
present_future_compare <- 100*nrow(present_future_compare)/nrow(Marxan_present)

df <- c(100, present_future_compare)

jpeg(paste0("./Marxan_future/Sensitivity/Present_Future/", Depth_class, "/SSOLN overlap Present_Future top ", Top_percent*100, " percent.jpg"))
barplot(df, xlab = c("Time period"), ylab = paste0("Percent similar between time periods for pld ", Baseline_pld))
dev.off()

paste0("./Marxan_future/Sensitivity/", Depth_class, "/Present_Future/SSOLN overlap Present_Future top ", Top_percent*100, " percent.jpg")
#####
####
###
##
#END
