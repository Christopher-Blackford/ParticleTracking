##############################################################################################################################
##############################################################################################################################
###Overlapping depth classes

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
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_study_extent_you_will_be_using.R")

Intertidal <- readOGR("./Marxan/To_shapefile/Intertidal/pld22", "MarxIntertidal_pld22")

Nearshore <- readOGR("./Marxan/To_shapefile/Nearshore/pld56", "MarxNearshore_pld56")

Offshore <- readOGR("./Marxan/To_shapefile/Offshore/pld48", "MarxOffshore_pld48")

All_marxan <- merge(Intertidal@data, Nearshore@data, by = "Poly_ID")
All_marxan <- merge(All_marxan, Offshore@data, by = "Poly_ID")

All_marxan <- dplyr::rename(All_marxan, Int_SSOLN = SSOLN.x, Near_SSOLN = SSOLN.y, Off_SSOLN = SSOLN)
All_marxan[is.na(All_marxan)]<-0

All_marxan$SUMM_SSOLN <- All_marxan$Int_SSOLN + All_marxan$Near_SSOLN + All_marxan$Off_SSOLN

#Getting columns to describe if they are in multiple depth layers
All_marxan$Int_Near <- All_marxan$Int_SSOLN * All_marxan$Near_SSOLN
All_marxan$Int_Off <- All_marxan$Int_SSOLN * All_marxan$Off_SSOLN
All_marxan$Near_Off <- All_marxan$Near_SSOLN * All_marxan$Off_SSOLN
  
All_SSOLN <- merge(ConPoly, All_marxan)



dir.create("./Marxan/To_shapefile/All_overlap")

writeOGR(All_SSOLN, dsn = "./Marxan/To_shapefile/All_overlap", layer = paste0("Marxan_all_depth_class"),
                driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

###############
#Getting top 10%
Top_percent <- 0.1

Top_percent_df <- All_marxan
Top_percent_df$included <- 1
Top_percent_df$RAND <- sample(1:nrow(Top_percent_df), nrow(Top_percent_df), replace = FALSE)
Top_percent_df <- Top_percent_df[with(Top_percent_df, order(-SUMM_SSOLN, RAND)),]
Top_percent_df <- Top_percent_df[c(1:ceiling(Top_percent*nrow(Top_percent_df))),]

Top_percent <- merge(ConPoly, Top_percent_df)

writeOGR(Top_percent, dsn = "K:/Christopher_PhD/temp", layer = paste0("All_depth_top_10"),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

