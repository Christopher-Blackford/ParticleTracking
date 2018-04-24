##############################################################################################################################
##############################################################################################################################
###R_output_to_Marxan.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

###READ.ME:
#This file takes: 
#Output files from Present_particle_locations.R of intertidal/nearshore/offshore data and
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
###[1] Load in Coastal BC shapefile - identical to - [2] Setting up study extent you will be using to clip your larval release points to your BC study extent in  "Present_Particle_locations"
###[2] Generate Marxan core file from Particle_tracking output
###[3] Preparing Marxan input files

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
Depth_class <- "Intertidal" #choices are "Intertidal", "Nearshore", "Offshore"
pld <- 22
my_resolution <- 10000 #defines raster cell size and controls for biased larvae release
target_percent_of_total <- 0.25

#Boundary file trickyness
Include_boundary_file <- TRUE
Edge_with_itself <- TRUE

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Load in Coastal BC shapefile - identical to - [2] Setting up study extent you will be using to clip your larval release points to your BC study extent in  "Present_Particle_locations"

source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_study_extent_you_will_be_using.R")
rm(My_BC_projection, NAD_projection, r)

########################################################################
########################################################################
########################################################################
########################################################################
###[2] Generate Marxan core file from Particle_tracking output
#Read in Particle_tracking output
R_file <- read.csv(paste0("./output_keep/Con_df/", Depth_class, "/pld", pld, "/", Depth_class, "_pld", pld, ".csv"))

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

input_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan/From_R/input/", Depth_class)
dir.create(input_directory)
input_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan/From_R/input/", Depth_class, "/pld", pld)
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
pu_file$cost <- 1

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

spf <- 1000 #"Species penalty factor" - penalty for not achieving a target
spf <- rep(spf, times = length(target))

species_file <- data.frame(id = c(1,2),
                           name = name,
                           #type = 1, Ignore type since you aren't using block file anymore
                           target = target,
                           spf = spf)

write.table(species_file, file=paste0(input_directory, "/species.dat"), row.names=FALSE, sep=",", quote=FALSE)

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
#Boundary_length <- Boundary_length[c("id1", "id2", "Freq")]


#Adding boundaries on edges
if (Edge_with_itself == TRUE){
  
  for (repeat_8 in 1:8){ #FIX THIS
  for (i in 1:nrow(Boundary_length)){
  if (nrow(Boundary_length[Boundary_length$Min_cell == i ,]) + 
      nrow(Boundary_length[Boundary_length$Max_cell == i ,]) < 8){
    Boundary_length <- rbind(i, Boundary_length)}
  }
    print(repeat_8)
}
 
 Boundary_length$Freq <- 1
  Boundary_length <- Boundary_length[c("id1", "id2", "Freq")]
  Boundary_length <- Boundary_length[with(Boundary_length, order(id1, id2)), ]
  row.names(Boundary_length) <- 1:nrow(Boundary_length)
}

write.table(Boundary_length, file=paste0(input_directory, "/boundary.dat"), row.names=FALSE, sep=",", quote=FALSE)

rm(i)
}
########################################################################
########################################################################
########################################################################
#Block definition file
#ignore block file for now - testing out summed conservation target feature approach
#block_file <- data.frame(type = 1, prop = proportion_protected, spf = 10)
#write.table(block_file, file=paste0(directory, "block.dat"), row.names=FALSE, sep=",", quote=FALSE)


########################################################################
########################################################################
########################################################################
########################################################################
###Setting up Marxan Parameters
input_directions <- read.delim("./Marxan/From_R/input.dat", header = TRUE, as.is = TRUE)

input_directory_marxan <- gsub(x = input_directory, pattern = "/", replacement = "\\\\", fixed = TRUE)

output_directory <- paste0("K:/Christopher_PhD/Github/ParticleTracking/Marxan/From_R/output/", Depth_class)
dir.create(output_directory)
output_directory <- paste0(output_directory, "/pld", pld)
dir.create(output_directory)
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
#Getting Marxan shapefile output

Marxan_output <- read.table(paste0(output_directory, "/mar_out_ssoln.dat"))
Marxan_output <- rename(Marxan_output, c("V1" = "Poly_ID", "V2" = "SSOLN"))

Marxan_spatial <- sp::merge(ConPoly, Marxan_output, by = "Poly_ID", all.x = TRUE)

shapefile_directory <- paste0("./Marxan/To_shapefile/", Depth_class, "/pld", pld)
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

