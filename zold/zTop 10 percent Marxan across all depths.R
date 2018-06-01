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

All_depth <- read.table("./Marxan/From_R/output/hexagon/All/Mar_out_ssoln.dat")
All_depth <- dplyr::rename(All_depth, Poly_ID = V1, SSOLN = V2)

All_depth <- merge(All_depth, ConPoly, by = "Poly_ID")

#Taking top percent
Marxan_output <- All_depth
Marxan_output$included <- 1
Marxan_output$RAND <- sample(1:nrow(Marxan_output), nrow(Marxan_output), replace = FALSE)
Marxan_output <- Marxan_output[with(Marxan_output, order(-SSOLN, RAND)),]
Marxan_output <- Marxan_output[c(1:ceiling(Top_percent*nrow(Marxan_output))),]

All_Top_10 <- sp::merge(ConPoly, Marxan_output, by = "Poly_ID", all.x=FALSE)

writeOGR(All_Top_10, dsn = "./Marxan/Top_10percent/All", layer = "All_Top10", 
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)



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

