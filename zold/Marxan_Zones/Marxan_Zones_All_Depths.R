##############################################################################################################################
##############################################################################################################################
###Marxan_Zones_All_Depths.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

#
##
###
####
#####
#Clear workspace
rm(list=ls())

###################TABLE OF CONTENTS

###################Loading required packages:
require(data.table)
require(tidyverse)
require(rgdal)
require(rgeos)
require(maptools)
require(plyr)
require(spatialEco)

setwd("K:/Christopher_PhD/Github/ParticleTracking")

###################Loading functions:
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/functions/my_point_in_poly.R")

########################################################################
########################################################################
#Setting depth class and pld
Intertidal_pld <- 22
Nearshore_pld <- 56
Offshore_pld <- 48

Intertidal <- read.csv(paste0("./output_keep/Con_df/hexagon/Intertidal/pld", Intertidal_pld, "/Intertidal_pld", Intertidal_pld, ".csv"))
Nearshore <- read.csv(paste0("./output_keep/Con_df/hexagon/Nearshore/pld", Nearshore_pld, "/Nearshore_pld", Nearshore_pld, ".csv"))
Offshore <- read.csv(paste0("./output_keep/Con_df/hexagon/Offshore/pld", Offshore_pld, "/Offshore_pld", Offshore_pld, ".csv"))

my_resolution <- 10000 #defines raster cell size and controls for biased larvae release
target_percent_of_total <- 0.25

########################################################################
########################################################################
########################################################################
###[1] Load in Coastal BC shapefile - identical to - [2] Setting up study extent you will be using to clip your larval release points to your BC study extent in  "Present_Particle_locations"

source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")
ConPoly_df <- ConPoly@data
Ecozones <- readOGR("./cuke_present/StudyExtent/BC_EcozonesMask", "Ecozones_BC")
Ecozones_df <- Ecozones@data

Intertidal <- merge(ConPoly_df, Intertidal, by = "Poly_ID", all.x = TRUE)
Nearshore <- merge(ConPoly_df, Nearshore, by = "Poly_ID", all.x = TRUE)
Offshore <- merge(ConPoly_df, Offshore, by = "Poly_ID", all.x = TRUE)

rm(My_BC_projection, NAD_projection)

########################################################################
########################################################################
########################################################################
###[3] Preparing Marxan input files

input_directory <- paste0("./Marxan_Zones/From_R/input")

########################################################################
########################################################################
########################################################################
#Planning unit file
#Strangely you don't need to include spatial point location of where the planning units are??
pu_file <- ConPoly_df
pu_file$area <- 1
pu_file <- dplyr::rename(pu_file, id = Poly_ID)

write.table(pu_file, file=paste0(input_directory, "/pu.dat"), row.names=FALSE, sep=",", quote=FALSE)

########################################################################
########################################################################
########################################################################
#Feature file

Depth_class <- list(Intertidal, Nearshore, Offshore)
target <- NULL

for (i in 1:length(Depth_class)){
  sum_conservation <- Depth_class[[i]][c("mean_Local_retention", "mean_Eigenvector_centrality")]
  sum_conservation[is.na(sum_conservation)] <- 0
  
  sum_Local_retention <- dplyr::summarise(sum_conservation, sum(mean_Local_retention))
  sum_Eigenvector_centrality <- dplyr::summarise(sum_conservation, sum(mean_Eigenvector_centrality))
  
  temp <- as.numeric(c(target_percent_of_total*sum_Local_retention, target_percent_of_total*sum_Eigenvector_centrality))
  target <- append(target, temp)
}

name = c("Intertidal_LR_target", "Intertidal_EC_target", "Nearshore_LR_target", "Nearshore_EC_target", "Offshore_LR_target", "Offshore_EC_target")
fpf <- 1 #"Feature penalty factor" - penalty for not achieving a target
feature_file <- data.frame(id = 1:6, target = target, fpf = fpf,  name = name)

write.table(feature_file, file=paste0(input_directory, "/feat.dat"), row.names=FALSE, sep=",", quote=FALSE)
rm(temp, target, fpf, name, sum_conservation, sum_Local_retention, sum_Eigenvector_centrality)

########################################################################
########################################################################
########################################################################
#Planning unit vs feature file
Depth_class <- list(Intertidal, Nearshore, Offshore)
feature_list_id <- c(1,3,5)
puvfeat <- NULL

for (i in 1:length(Depth_class)){
  pu_v_conservation <- Depth_class[[i]][c("Poly_ID", "mean_Local_retention", "mean_Eigenvector_centrality")]
  pu_v_conservation[is.na(pu_v_conservation)] <- 0
  
  #For local retention
  puvfeat1 <- pu_v_conservation[c("Poly_ID", "mean_Local_retention")]
  puvfeat1$featureid <- feature_list_id[i]
  puvfeat1 <- dplyr::rename(puvfeat1, puid = Poly_ID, amount = mean_Local_retention)
  #For Eigenvector centrality
  puvfeat2 <- pu_v_conservation[c("Poly_ID", "mean_Eigenvector_centrality")]
  puvfeat2$featureid <- feature_list_id[i]+1
  puvfeat2 <- dplyr::rename(puvfeat2, puid = Poly_ID, amount = mean_Eigenvector_centrality)
  
  pu_v_conservation <- rbind(puvfeat1, puvfeat2)
  puvfeat <- rbind(puvfeat, pu_v_conservation)
}

puvfeat <- puvfeat[c("featureid", "puid", "amount")]

write.table(puvfeat, file=paste0(input_directory, "/puvfeat.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(feature_list_id, pu_v_conservation, puvfeat1, puvfeat2)

########################################################################
########################################################################
########################################################################
#Zones file
zones_file <- Ecozones_df
zones_file <- dplyr::rename(zones_file, zoneid = FID_CA_EZ_, zonename = Name)
zones_file <- zones_file[c("zoneid", "zonename")]
row.names(zones_file) <- 1:nrow(zones_file)
zones_file$zonename <- gsub(x = zones_file$zonename, pattern = " ", replacement = "_", fixed = TRUE)

write.table(zones_file, file=paste0(input_directory, "/zones.dat"), row.names=FALSE, sep=",", quote=FALSE)


########################################################################
########################################################################
########################################################################
#Costs file
costs <- data.frame(costid = 1, costname = "area")

write.table(costs, file=paste0(input_directory, "/costs.dat"), row.names=FALSE, sep=",", quote=FALSE)

########################################################################
########################################################################
########################################################################
#Zone costs file
zonecost <- data.frame(zoneid = zones_file$zoneid, costid = costs$costid, multiplier = 1)

write.table(zonecost, file=paste0(input_directory, "/zonecost.dat"), row.names=FALSE, sep=",", quote=FALSE)

########################################################################
########################################################################
########################################################################
#Boundary length file

Boundary_length <- gTouches(ConPoly, byid = TRUE) #tells you if a cell is touching another cell
Boundary_length <- as.table(Boundary_length)
Boundary_length <- as.data.frame(Boundary_length)

Boundary_length <- Boundary_length[Boundary_length$Freq == "TRUE", ]
Boundary_length$Freq <- 1
Boundary_length <- dplyr::rename(Boundary_length, id1 = Var1, id2 = Var2)
Boundary_length$Min_cell <- NA ; Boundary_length$Max_cell <- NA

Boundary_length$id1 <- as.numeric(Boundary_length$id1) ; Boundary_length$id2 <- as.numeric(Boundary_length$id2)
Boundary_length <- Boundary_length[with(Boundary_length, order(id1, id2)), ]
row.names(Boundary_length) <- 1:nrow(Boundary_length)


for (i in 1:nrow(Boundary_length)){
  Boundary_length$Min_cell[i] <- as.character(min(c(Boundary_length$id1[i], Boundary_length$id2[i])))
  Boundary_length$Max_cell[i] <- as.character(max(c(Boundary_length$id1[i], Boundary_length$id2[i])))
}

Boundary_length$Unique_cell <- paste0(Boundary_length$Min_cell, Boundary_length$Max_cell)
Boundary_length <- Boundary_length[!duplicated(Boundary_length$Unique_cell) ,]
row.names(Boundary_length) <- 1:nrow(Boundary_length)
Boundary_length <- Boundary_length[c("id1", "id2", "Freq")]

Polygon_list <- unique(c(Boundary_length$id1, Boundary_length$id2))

#Adding boundaries on edges
for (i in Polygon_list){
  if (nrow(Boundary_length[Boundary_length$id1 == i ,]) + nrow(Boundary_length[Boundary_length$id2 == i ,]) < 6){
    
    num_connections <- nrow(Boundary_length[Boundary_length$id1 == i ,]) + nrow(Boundary_length[Boundary_length$id2 == i ,])
    for(repeat_number in 1:(6-num_connections)){Boundary_length <- rbind(i, Boundary_length)}
  }
}
Boundary_length$Freq <- 1
Boundary_length <- dplyr::rename(Boundary_length, boundary = Freq)
Boundary_length <- Boundary_length[with(Boundary_length, order(id1, id2)), ]
row.names(Boundary_length) <- 1:nrow(Boundary_length)


write.table(Boundary_length, file=paste0(input_directory, "/bounds.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(i)


########################################################################
########################################################################
########################################################################
#Planning unit lock file - sets pu's to a single zone

#Get centroid points of each cell since Marxan requires point data
Centroid_points <- SpatialPointsDataFrame(gCentroid(ConPoly, byid = TRUE), ConPoly@data, match.ID = "Poly_ID")
Centroid_points <- spTransform(Centroid_points, ConPoly@proj4string)

pulock_file <- data.frame(sp::over(Centroid_points, Ecozones))
pulock_file <- pulock_file[c("FID_CA_EZ_", "Name")]
pulock_file$Poly_ID <- row.names(pulock_file)
pulock_file <- dplyr::rename(pulock_file, puid = Poly_ID, zoneid = FID_CA_EZ_)
pulock_file <- pulock_file[c("puid", "zoneid")]

write.table(pulock_file, file=paste0(input_directory, "/pulock.dat"), row.names=FALSE, sep=",", quote=FALSE)



########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
#Running Marxan

system(paste("K:/Christopher_PhD/Github/ParticleTracking/Marxan_Zones/From_R/MarZone_x64.exe K:/Christopher_PhD/Github/ParticleTracking/Marxan_Zones/From_R/input.dat", sep = ""))


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################



#####
####
###
##
#END