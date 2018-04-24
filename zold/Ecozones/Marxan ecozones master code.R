##############################################################################################################################
##############################################################################################################################
###Marxan ecozones master code

#
##
###
####
#####
#Clear workspace
rm(list=ls())

###################Loading required packages:
require(data.table)
require(tidyverse)
require(rgdal)
require(rgeos)
require(maptools)
require(plyr)
require(spatialEco)

setwd("K:/Christopher_PhD/Github/ParticleTracking")

########################################################################
########################################################################
#Setting depth class and pld

Ecozone_region <- "Northern_Shelf"

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

###Clipping to spatial extent
if (Ecozone_region == "Southern_Shelf"){
  Ecozones <- Ecozones[Ecozones$FID_CA_EZ_ == 1,] #Southern_Shelf
  } else if (Ecozone_region == "Straight_of_Georgia"){
  Ecozones <- Ecozones[Ecozones$FID_CA_EZ_ == 2,] #Straight_of_Georgia
  } else if (Ecozone_region == "Northern_Shelf"){
  Ecozones <- Ecozones[Ecozones$FID_CA_EZ_ == 6,] #Northern_Shelf
  }

ConPoly <- ConPoly[Ecozones,]
ConPoly_df <- ConPoly@data

Intertidal <- merge(ConPoly_df, Intertidal, by = "Poly_ID", all.x = TRUE)
Nearshore <- merge(ConPoly_df, Nearshore, by = "Poly_ID", all.x = TRUE)
Offshore <- merge(ConPoly_df, Offshore, by = "Poly_ID", all.x = TRUE)

rm(My_BC_projection, NAD_projection)

########################################################################
########################################################################
########################################################################
###[2] Generate Marxan core file from Particle_tracking output

#Get centroid points of each cell since Marxan requires point data
Centroid_points <- SpatialPointsDataFrame(gCentroid(ConPoly, byid = TRUE), ConPoly@data, match.ID = "Poly_ID")
Centroid_points <- data.frame(Centroid_points@coords)
Centroid_points$Poly_ID <- row.names(Centroid_points)
#Merge with ConPoly to get point data combined with local retention and centrality
Study_extent <- sp::merge(ConPoly, Centroid_points, by = "Poly_ID")

Study_extent_file <- Study_extent@data
rm(Centroid_points)

########################################################################
########################################################################
########################################################################
###[3] Preparing Marxan input files

input_directory <- paste0("./Marxan/From_R/input/hexagon/Ecozones/", Ecozone_region)

########################################################################
########################################################################
########################################################################
#Planning unit file
pu_file <- Study_extent_file

pu_file <- rename(pu_file, c("Poly_ID" = "id","x" = "xloc", "y" = "yloc"))
pu_file$cost <- 1 #cost is the same for all the cells
pu_file$status <- 0 #any cell could make up reserve solution
pu_file <- pu_file[c("id", "cost", "status", "xloc", "yloc")]

write.table(pu_file, file=paste0(input_directory, "/pu.dat"), row.names=FALSE, sep=",", quote=FALSE)



########################################################################
########################################################################
#Conservation feature file (species file)
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
spf <- 100 #"Species penalty factor" - penalty for not achieving a target
species_file <- data.frame(id = 1:6, name = name, target = target, spf = spf)

write.table(species_file, file=paste0(input_directory, "/species.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(temp, name, sum_conservation, sum_Local_retention, sum_Eigenvector_centrality)


########################################################################
########################################################################
########################################################################
#Planning unit versus Conservation feature file (puvspr2 file)
Depth_class <- list(Intertidal, Nearshore, Offshore)
species_list_id <- c(1,3,5)
puvspr2_file <- NULL

for (i in 1:length(Depth_class)){
  pu_v_conservation <- Depth_class[[i]][c("Poly_ID", "mean_Local_retention", "mean_Eigenvector_centrality")]
  pu_v_conservation[is.na(pu_v_conservation)] <- 0
  
  #For local retention
  puvspr2_file1 <- pu_v_conservation[c("Poly_ID", "mean_Local_retention")]
  puvspr2_file1$species <- species_list_id[i]
  puvspr2_file1 <- dplyr::rename(puvspr2_file1, pu = Poly_ID, amount = mean_Local_retention)
  #For Eigenvector centrality
  puvspr2_file2 <- pu_v_conservation[c("Poly_ID", "mean_Eigenvector_centrality")]
  puvspr2_file2$species <- species_list_id[i]+1
  puvspr2_file2 <- dplyr::rename(puvspr2_file2, pu = Poly_ID, amount = mean_Eigenvector_centrality)
  
  pu_v_conservation <- rbind(puvspr2_file1, puvspr2_file2)
  puvspr2_file <- rbind(puvspr2_file, pu_v_conservation)
}

write.table(puvspr2_file, file=paste0(input_directory, "/puvspr2.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(pu_v_conservation, puvspr2_file1, puvspr2_file2)



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
Boundary_length <- Boundary_length[with(Boundary_length, order(id1, id2)), ]
row.names(Boundary_length) <- 1:nrow(Boundary_length)


write.table(Boundary_length, file=paste0(input_directory, "/boundary.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(i)


########################################################################
########################################################################
###Setting up Marxan Parameters
input_directions <- read.delim("./Marxan/From_R/input.dat", header = TRUE, as.is = TRUE)

input_directory_marxan <- gsub(x = input_directory, pattern = "/", replacement = "\\\\", fixed = TRUE)

output_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan/From_R/output/hexagon/Ecozones/", Ecozone_region)
output_directory_marxan <- gsub(x = output_directory, pattern = "/", replacement = "\\\\", fixed = TRUE)

input_edit <- data.frame(lapply(input_directions, function(x) {gsub(pattern = "INPUTDIR.*", replacement = "INPUTDIR", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "INPUTDIR.*", replacement = paste0("INPUTDIR ", input_directory_marxan),  x = x)}))

input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "OUTPUTDIR.*", replacement = "OUTPUTDIR", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "OUTPUTDIR.*", replacement = paste0("OUTPUTDIR ", output_directory_marxan), x = x)}))

write.table(input_edit, file= "./Marxan/From_R/input.dat", row.names=FALSE, sep=",", quote=FALSE)

rm(input_directions, input_directory_marxan)


########################################################################
########################################################################
########################################################################
########################################################################
#Running Marxan

system(paste("K:/Christopher_PhD/Github/ParticleTracking/Marxan/From_R/Marxan_x64.exe K:/Christopher_PhD/Github/ParticleTracking/Marxan/From_R/input.dat", sep = ""))

########################################################################
########################################################################
########################################################################
########################################################################






########################################################################
########################################################################
########################################################################
#Getting Marxan shapefile output

Marxan_output <- read.table(paste0(output_directory, "/mar_out_ssoln.dat"))
Marxan_output <- rename(Marxan_output, c("V1" = "Poly_ID", "V2" = "SSOLN"))

Marxan_spatial <- sp::merge(ConPoly, Marxan_output, by = "Poly_ID", all.x = TRUE)

shapefile_directory <- paste0("./Marxan/To_shapefile/hexagon/Ecozone/", Ecozone_region)

writeOGR(Marxan_spatial, dsn = shapefile_directory, layer = paste0("Int", Intertidal_pld, "_Near", Nearshore_pld, "_Off", Offshore_pld),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

#####
####
###
##
#END









