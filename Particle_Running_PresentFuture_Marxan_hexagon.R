##############################################################################################################################
##############################################################################################################################
###R_output_to_Marxan.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

###READ.ME:
#This file takes: 
#Output files from Future_particle_hexagon_locations.R of intertidal/nearshore/offshore data and
#transforms it into Marxan-ready files

#
##
###
####
#####
#Clear workspace
rm(list=ls())
marxan.time <- proc.time() # for one run-through

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

########################################################################
########################################################################
#Setting depth class and pld

#Intertidal PLDs: 6, 9, 21, 24, 28, 29, 40                 #+/- 2 days around average 22: 22, 20, 23
#Nearshore PLDs: 3, 4, 6, 9, 24, 31, 38, 45, 52, 58, 60, 61, 78, 91, 95, 105, 109, 120                  #+/- 2 days around average 56: 56, 54, 55, 57
#Offshore PLDs: 27, 45, 90                                #+/- 2 days around average 48: 48, 46, 47, 49, 50

Depth_class <- "Intertidal" #choices are "Intertidal", "Nearshore", "Offshore"
pld <- 2
my_resolution <- 10000 #defines raster cell size and controls for biased larvae release
target_percent_of_total <- 0.25

cell_cost <- 0.1
boundary_length_modifier <- 1
boundary_length_cost <- 0.01
spf <- 100
cost_threshold <- 12
species_missing_threshold <- 1

#Boundary file trickyness
Include_boundary_file <- TRUE
Edge_with_itself <- TRUE

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Load in Coastal BC shapefile - identical to - [2] Setting up study extent you will be using to clip your larval release points to your BC study extent in  "Present_Particle_locations"

source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")
rm(My_BC_projection, NAD_projection)

########################################################################
########################################################################
########################################################################
########################################################################
###[2] Generate Marxan core file from Particle_tracking output
#Read in Particle_tracking output
R_file <- read.csv(paste0("./output_keep_future/Con_df/hexagon/", Depth_class, "/pld", pld, "/", Depth_class, "_pld", pld, ".csv"))

#Get centroid points of each cell since Marxan requires point data
Centroid_points <- SpatialPointsDataFrame(gCentroid(ConPoly, byid = TRUE), ConPoly@data, match.ID = "Poly_ID")
Centroid_points <- data.frame(Centroid_points@coords)
Centroid_points$Poly_ID <- row.names(Centroid_points)
#Merge with ConPoly to get point data combined with local retention and centrality
Marxan_spatial <- sp::merge(ConPoly, Centroid_points, by = "Poly_ID")
Marxan_spatial <- sp::merge(Marxan_spatial, R_file, by = "Poly_ID", all.x = TRUE)

Marxan_file <- Marxan_spatial@data
rm(Centroid_points, R_file)
########################################################################
########################################################################
########################################################################
########################################################################
###[3] Preparing Marxan input files

input_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan_future/From_R/input/hexagon/", Depth_class)
#dir.create(input_directory)
input_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan_future/From_R/input/hexagon/", Depth_class, "/pld", pld)
dir.create(input_directory)

########################################################################
########################################################################
########################################################################
#Planning unit file
pu_file <- Marxan_file

###Ensuring areas outside Depth Class aren't included in reserve
for (i in 1:nrow(pu_file)){
  if (is.na(pu_file$mean_Local_retention[i]) == TRUE) {
    pu_file$status[i] <- 3
  }
  else {pu_file$status[i] <- 0}
} 
#####

pu_file <- rename(pu_file, c("Poly_ID" = "id","x" = "xloc", "y" = "yloc"))
pu_file$cost <- cell_cost
pu_file <- pu_file[c("id", "cost", "status", "xloc", "yloc")]

write.table(pu_file, file=paste0(input_directory, "/pu.dat"), row.names=FALSE, sep=",", quote=FALSE)


########################################################################
########################################################################
########################################################################
#Conservation feature file (species file)
sum_conservation <- Marxan_file[c("mean_Local_retention", "mean_Eigenvector_centrality")]
sum_conservation[is.na(sum_conservation)] <- 0

sum_Local_retention <- dplyr::summarise(sum_conservation, sum(mean_Local_retention))
sum_Eigenvector_centrality <- dplyr::summarise(sum_conservation, sum(mean_Eigenvector_centrality))

target <- as.numeric(c(target_percent_of_total*sum_Local_retention, target_percent_of_total*sum_Eigenvector_centrality))
name = c("Local_retention_target", "Eigenvector_centrality_target")

spf <- 100 #"Species penalty factor" - penalty for not achieving a target
spf <- rep(spf, times = length(target))

species_file <- data.frame(id = c(1,2),
                           name = name,
                           target = target,
                           spf = spf)

species_file <- species_file[(c("id", "target", "spf", "name"))]

write.table(species_file, file=paste0(input_directory, "/spec.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(name, sum_conservation)
########################################################################
########################################################################
########################################################################
#Planning unit versus Conservation feature file (puvspr2 file)
puvspr2_file <- Marxan_file[c("Poly_ID", "mean_Local_retention", "mean_Eigenvector_centrality")]
puvspr2_file[is.na(puvspr2_file)] <- 0

#For local retention
puvspr2_file1 <- puvspr2_file[c("Poly_ID", "mean_Local_retention")]
puvspr2_file1$species <- 1
puvspr2_file1 <- dplyr::rename(puvspr2_file1, pu = Poly_ID, amount = mean_Local_retention)
#For Eigenvector centrality
puvspr2_file2 <- puvspr2_file[c("Poly_ID", "mean_Eigenvector_centrality")]
puvspr2_file2$species <- 2
puvspr2_file2 <- dplyr::rename(puvspr2_file2, pu = Poly_ID, amount = mean_Eigenvector_centrality)

puvspr2_file <- rbind(puvspr2_file1, puvspr2_file2)
puvspr2_file <- puvspr2_file[c("species", "pu", "amount")]
puvspr2_file <- arrange(puvspr2_file, -desc(pu))

write.table(puvspr2_file, file=paste0(input_directory, "/puvspr2.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(puvspr2_file1, puvspr2_file2)
########################################################################
########################################################################
########################################################################
#Boundary length file
if (Include_boundary_file == TRUE){

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
if (Edge_with_itself == TRUE){
  
for (i in Polygon_list){
  if (nrow(Boundary_length[Boundary_length$id1 == i ,]) + nrow(Boundary_length[Boundary_length$id2 == i ,]) < 6){
    
    num_connections <- nrow(Boundary_length[Boundary_length$id1 == i ,]) + nrow(Boundary_length[Boundary_length$id2 == i ,])
    for(repeat_number in 1:(6-num_connections)){Boundary_length <- rbind(i, Boundary_length)}
    }
}
  Boundary_length$boundary <- boundary_length_cost
  Boundary_length <- Boundary_length[with(Boundary_length, order(id1, id2)), ]
  row.names(Boundary_length) <- 1:nrow(Boundary_length)
  Boundary_length <- Boundary_length[c("id1", "id2", "boundary")]
  
}

write.table(Boundary_length, file=paste0(input_directory, "/boundary.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(i)
}


########################################################################
########################################################################
########################################################################
########################################################################
###Setting up Marxan Parameters
input_directions <- read.delim("./Marxan_future/From_R/input.dat", header = TRUE, as.is = TRUE)

input_directory_marxan <- gsub(x = input_directory, pattern = "/", replacement = "\\\\", fixed = TRUE)

output_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan_future/From_R/output/hexagon/", Depth_class)
dir.create(output_directory)
output_directory <- paste0(output_directory, "/pld", pld)
dir.create(output_directory)
output_directory_marxan <- gsub(x = output_directory, pattern = "/", replacement = "\\\\", fixed = TRUE)

#Modifying input directory
input_edit <- data.frame(lapply(input_directions, function(x) {gsub(pattern = "INPUTDIR.*", replacement = "INPUTDIR", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "INPUTDIR.*", replacement = paste0("INPUTDIR ", input_directory_marxan),  x = x)}))
#Modifying output directory
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "OUTPUTDIR.*", replacement = "OUTPUTDIR", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "OUTPUTDIR.*", replacement = paste0("OUTPUTDIR ", output_directory_marxan), x = x)}))
#Modifying boundary length modifier
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "BLM.*", replacement = "BLM", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "BLM*", replacement = paste0("BLM ", boundary_length_modifier), x = x)}))
#Modifying threshold
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "COSTTHRESH.*", replacement = "COSTTHRESH", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "COSTTHRESH.*", replacement = paste0("COSTTHRESH ", cost_threshold), x = x)}))
#Modifying species missing if proportion is under "x"
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "MISSLEVEL.*", replacement = "MISSLEVEL.", x = x)}))
input_edit <- data.frame(lapply(input_edit, function(x) {gsub(pattern = "MISSLEVEL.*", replacement = paste0("MISSLEVEL ", species_missing_threshold), x = x)}))

write.table(input_edit, file= "./Marxan_future/From_R/input.dat", row.names=FALSE, sep=",", quote=FALSE)

rm(input_directions, input_directory_marxan)
########################################################################
########################################################################
########################################################################
########################################################################
#Running Marxan

system(paste("K:/Christopher_PhD/Github/ParticleTracking/Marxan_future/From_R/Marxan_x64.exe K:/Christopher_PhD/Github/ParticleTracking/Marxan_future/From_R/input.dat", sep = ""))

########################################################################
########################################################################
########################################################################
########################################################################
#Getting Marxan shapefile output

Marxan_output <- read.table(paste0(output_directory, "/mar_out_ssoln.dat"))
Marxan_output <- rename(Marxan_output, c("V1" = "Poly_ID", "V2" = "SSOLN"))

Marxan_spatial <- sp::merge(ConPoly, Marxan_output, by = "Poly_ID", all.x = TRUE)

shapefile_directory <- paste0("./Marxan_future/To_shapefile/hexagon/", Depth_class, "/pld", pld)
dir.create(shapefile_directory)

writeOGR(Marxan_spatial, dsn = shapefile_directory, layer = paste0("Marx", Depth_class, "_pld", pld),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)


proc.time() - marxan.time
##########
#########
########
#######
######
#####
####
###
##
#END

Marxan_sum <- read.table(paste0(output_directory, "/mar_out_sum.dat"))
