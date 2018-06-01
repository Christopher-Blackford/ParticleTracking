##############################################################################################################################
##############################################################################################################################
###Top 10 percent for each depth class

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
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")

Top_percent <- 0.1

Intertidal <- read.table("./Marxan/From_R/output/hexagon/Intertidal/pld22/Mar_out_ssoln.dat")
Nearshore <- read.table("./Marxan/From_R/output/hexagon/Nearshore/pld56/Mar_out_ssoln.dat")
Offshore <- read.table("./Marxan/From_R/output/hexagon/Offshore/pld48/Mar_out_ssoln.dat")

Habitat_classes <- list(Intertidal, Nearshore, Offshore)
Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")

for(i in 1:length(Habitat_classes)){
  temp <- dplyr::rename(Habitat_classes[[i]], Poly_ID = V1, SSOLN = V2)
  assign(Habitat_classes_names[i], temp)}; rm(temp)

#Merging with actual Depth class extent
Intertidal_extent <- readOGR("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/Intertidal", "Intertidal")
Nearshore_extent <- readOGR("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/Nearshore", "Nearshore")
Offshore_extent <- readOGR("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/Offshore", "Offshore")

Intertidal <- merge(Intertidal, Intertidal_extent, by = "Poly_ID")
Nearshore <- merge(Nearshore, Nearshore_extent, by = "Poly_ID")
Offshore <- merge(Offshore, Offshore_extent, by = "Poly_ID")

Habitat_classes <- list(Intertidal, Nearshore, Offshore)
Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")

#Gets you dataframes showing top 10% for each depth class
for (i in 1:length(Habitat_classes_names)){
  
  Marxan_output <- Habitat_classes[[i]]
  #Marxan_output <- dplyr::rename(Marxan_output, Poly_ID = V1, SSOLN = V2)
  Marxan_output$included <- 1
  Marxan_output$RAND <- sample(1:nrow(Marxan_output), nrow(Marxan_output), replace = FALSE)
  Marxan_output <- Marxan_output[with(Marxan_output, order(-SSOLN, RAND)),]
  Marxan_output <- Marxan_output[c(1:ceiling(Top_percent*nrow(Marxan_output))),]
  
  assign(Habitat_classes_names[i], Marxan_output)
  }

Intertidal_top10 <- sp::merge(Intertidal_extent, Intertidal, by = "Poly_ID", all.x=FALSE)
Nearshore_top10 <- sp::merge(Nearshore_extent, Nearshore, by = "Poly_ID", all.x=FALSE)
Offshore_top10 <- sp::merge(Offshore_extent, Offshore, by = "Poly_ID", all.x=FALSE)

Habitat_classes_top10 <- list(Intertidal_top10, Nearshore_top10, Offshore_top10)

for (i in 1:length(Habitat_classes_top10)){
writeOGR(Habitat_classes_top10[[i]], dsn = paste0("./Marxan/Top_10percent/", Habitat_classes_names[i]), layer = paste0(Habitat_classes_names[i], "_top10"),
        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)}


#####
####
###
##
#END







####Overlaps between depth classes
#Int_Near <- merge(Intertidal, Nearshore, by = "Poly_ID")
#Int_Near_sp <- sp::merge(ConPoly, Int_Near, by = "Poly_ID", all.x = TRUE)
#writeOGR(Int_Near_sp, dsn = "./Marxan/To_shapefile/All_overlap", layer = "Int_Near_top10",
 #        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)


#Int_Off <- merge(Intertidal, Offshore, by = "Poly_ID")
#Int_Off_sp <- sp::merge(ConPoly, Int_Off, by = "Poly_ID", all.x = TRUE)
#writeOGR(Int_Off_sp, dsn = "./Marxan/To_shapefile/All_overlap", layer = "Int_Off_top10",
 #        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)


#Near_Off <- merge(Nearshore, Offshore, by = "Poly_ID")
#Near_Off_sp <- sp::merge(ConPoly, Near_Off, by = "Poly_ID", all.x = TRUE)
#writeOGR(Near_Off_sp, dsn = "./Marxan/To_shapefile/All_overlap", layer = "Near_Off_top10",
 #        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

