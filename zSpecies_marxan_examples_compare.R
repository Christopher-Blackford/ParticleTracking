##############################################################################################################################
##############################################################################################################################
###Similarity between species

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

setwd("K:/Christopher_PhD/Github/ParticleTracking") 

###################Comparing any 2

#Species_1 <- readOGR("./Marxan/To_shapefile/hexagon/Nearshore/pld3")
#Species_2 <- readOGR("./Marxan/To_shapefile/hexagon/Nearshore/pld120")

#Species_1 <- Species_1@data
#Species_2 <- Species_2@data

#Species_1$RAND <- sample(1:nrow(Species_1), nrow(Species_1), replace = FALSE)
#Species_1 <- Species_1[with(Species_1, order(-SSOLN, RAND)),]
#Species_1 <- Species_1[c(1:round(0.1*nrow(Species_1))),]
#Species_1 <- Species_1[with(Species_1, order(Poly_ID)),]
#Species_1 <- Species_1[,1]


#Species_2$RAND <- sample(1:nrow(Species_2), nrow(Species_2), replace = FALSE)
#Species_2 <- Species_2[with(Species_2, order(-SSOLN, RAND)),]
#Species_2 <- Species_2[c(1:round(0.1*nrow(Species_2))),]
#Species_2 <- Species_2[with(Species_2, order(Poly_ID)),]
#Species_2 <- Species_2[,1]


#over_lap <- length(intersect(Species_1, Species_2))/length(Species_1)


###################Comparing Current to Future all depth classes
Species_1 <- readOGR("./Thesis_Figures/StudyRegion/Data/zCurrent_Future", "Current_25")
Species_2 <- readOGR("./Thesis_Figures/StudyRegion/Data/zCurrent_Future", "Future_25")

Species_1 <- Species_1@data
Species_2 <- Species_2@data

Species_1 <- Species_1[,1]
Species_2 <- Species_2[,1]

over_lap <- length(intersect(Species_1, Species_2))/length(Species_1)
