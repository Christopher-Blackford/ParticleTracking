#############################################################################################################################
##############################################################################################################################
###Connectivity_Stephen_Rockfish.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

#
##
###
####
#####


#Clear workspace
rm(list=ls())
full.run.time <- proc.time() # 33 minutes

###################Loading required packages:
require(plyr)
require(data.table)
require(tidyverse)
require(rgdal)
require(rgeos)
require(maptools)
require(raster)
setwd("K:/Christopher_PhD/Github/ParticleTracking")

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Setting up study extent you will be using to clip your larval release points to your BC study extent


########################################################################
########################################################################
#Clipping to your study extent
Ecozone_mask <- readOGR("./cuke_present/StudyExtent/BC_EcozonesMask", "Ecozone_mask")
My_BC_projection <- proj4string(Ecozone_mask)

#Loading Remi's grid where larvae were released
grid <- readOGR("./cuke_present/StudyExtent/Starting_grid", "grid")
NAD_projection <- proj4string(grid)

#Dissolve into one polygon since so you can change grid dimensions
grid <- spTransform(grid, Ecozone_mask@proj4string) #For some reason not "identical" to My_BC_projection, check later
grid <- gUnaryUnion(grid)

#Intersecting - don't know why this works and ConPoly2 <- grid[Ecozone_mask,] doesn't
ConPoly <- gIntersection(grid, Ecozone_mask, byid = FALSE, drop_lower_td = TRUE) #This works, but you'll have to choose a shapefile that includes islands and doesn't cut-off at rivers 

########################################################################
########################################################################
#Bathymetry
Bathymetry_raster <- raster("./Connectivity_between_MPA_Stephen/Rockfish_project/input_data/bcmca_eco_physical_bathymetry_data/raster/bathymetry")

m <- c(-200, Inf,  1, -Inf, -201, NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

Bathymetry <- reclassify(Bathymetry_raster, rclmat)

Bathymetry <- projectRaster(Bathymetry, crs = ConPoly@proj4string)
Bathymetry <- intersect(Bathymetry, ConPoly)
Bathymetry <- aggregate(Bathymetry, fact=20, fun=mean, na.rm=TRUE)

memory.limit(size=15000)
Bathymetry <- rasterToPolygons(Bathymetry, fun=NULL, n=4, na.rm=TRUE, digits=2, dissolve=TRUE)

rm(Bathymetry_raster, m, rclmat)
########################################################################
########################################################################
#Benthic data
Benthic <- readOGR("./Connectivity_between_MPA_Stephen/Rockfish_project/input_data/bcmca_eco_physical_benthicclasses_data", "BCMCA_ECO_Physical_BenthicClasses_DATA")
Benthic@data$merge <- c(1:nrow(Benthic@data))
df <- Benthic@data
df <- cbind(read.fwf(file = textConnection(as.character(df[,"habClass"])), 
                     widths = c(1, 1, 1, 1), colClasses = "character", 
                     col.names = c("DepthClass", "Col2", "Substrate", "Benthic_pos")))
df$merge <- c(1:nrow(df))

df$Substrate[df$Substrate == 1] <- "Mud"
df$Substrate[df$Substrate == 2] <- "Sand"
df$Substrate[df$Substrate == 3] <- "Hard"
df$Substrate[df$Substrate == 9] <- "Unknown"

df$Benthic_pos[df$Benthic_pos == 1] <- "Ridge"
df$Benthic_pos[df$Benthic_pos == 2] <- "Depression"
df$Benthic_pos[df$Benthic_pos == 3] <- "Flat"
df$Benthic_pos[df$Benthic_pos == 4] <- "Slope"

Benthic <- merge(Benthic, df, by = "merge")

#just to see
writeRaster(Benthic, filename= "./Connectivity_between_MPA_Stephen/Rockfish_project/output_rasters/habitat_quality/inputs/Benthic.tiff", proj4string = Benthic_ras@crs, overwrite=TRUE)

#Only keep landscapes that are hard
Benthic <- Benthic[which(Benthic@data$Substrate == "Hard"),]
head(Benthic@data)

for (i in 1:nrow(Benthic@data)){
  if (Benthic@data$Benthic_pos[i] == "Flat"){
    Benthic@data$Value[i] <- 1 
  }
  
  else if (Benthic@data$Benthic_pos[i] == "Depression"){
    Benthic@data$Value[i] <- 2 
  }
  
  else if (Benthic@data$Benthic_pos[i] == "Ridge" | Benthic@data$Benthic_pos[i] == "Slope"){
    Benthic@data$Value[i] <- 3 
  }  
}

Benthic_new <- spTransform(Benthic, Ecozone_mask@proj4string)
row.names(Benthic_new) <- as.character(Benthic_new@data$merge)

Benthic_new <- gBuffer(Benthic_new, byid=TRUE, width=0) #Need to do this to avoid ring intersection

dataframed <- proc.time() 
Benthic_intersect <- gIntersection(Benthic_new, Bathymetry, byid = TRUE) #40 minutes
proc.time() - dataframed

backup <- Benthic_intersect

df <- data.frame(row.names(backup))
df <- plyr::rename(df, c("row.names.backup." = "mergeID"))
df$rows <- row.names(df)


sf <- data.frame(do.call('rbind', strsplit(as.character(df$mergeID)," ",fixed=TRUE)))
sf$rows <- row.names(sf)
sf <- plyr::rename(sf, c("X1" = "merge", "X2" = "Col2"))

sf2 <- Benthic@data
sf <- merge(sf, sf2, by = "merge")


df <- merge(sf, df, by = "rows")
df <- df[c("mergeID", "merge", "Value")]

Benthic <- SpatialPolygonsDataFrame(backup, df, match.ID = "mergeID")
row.names(Benthic@data) <- 1:nrow(Benthic@data)


#Turn into raster
r <- raster(extent(Bathymetry))
projection(r) <- proj4string(Benthic)
res(r) <- 1000 #in metres, size of each cell (e.g. 5000 = 5000x5000m) CHANGE BACK TO 5000

Benthic_ras <- rasterize(Benthic, r, field = "Value") #5 seconds

#Need to replace NAs with zeros for raster calculator
Benthic_ras_calc <- Benthic_ras
Benthic_ras_calc[is.na(Benthic_ras_calc[])] <- 0 

writeRaster(Benthic_ras_calc, filename= "./Connectivity_between_MPA_Stephen/Rockfish_project/output_rasters/habitat_quality/inputs/Benthic.tiff", proj4string = Benthic_ras@crs, overwrite=TRUE)

rm(df, sf, sf2, backup)
########################################################################
########################################################################
#High Rugosity
High_Rugosity <- readOGR("./Connectivity_between_MPA_Stephen/Rockfish_project/input_data/bcmca_eco_physical_highrugosity_data", "BCMCA_ECO_Physical_HighRugosity_DATA")

#project
High_Rugosity <- spTransform(High_Rugosity, Ecozone_mask@proj4string)
High_Rugosity@data

#clip
High_Rugosity_data <- gIntersection(High_Rugosity, Bathymetry, byid = FALSE) #43 minutes


High_Rugosity_Value <- 1
High_Rugosity_Value <- as.data.frame(High_Rugosity_Value)

High_Rugosity <- SpatialPolygonsDataFrame(High_Rugosity_data, High_Rugosity_Value)
plot(High_Rugosity)

#Turn into raster
r <- raster(extent(Bathymetry))
projection(r) <- proj4string(High_Rugosity)
res(r) <- 1000 #in metres, size of each cell (e.g. 5000 = 5000x5000m) CHANGE BACK TO 5000

High_Rugosity_ras <- rasterize(High_Rugosity, r, field = "High_Rugosity_Value") #5 seconds

#Need to replace NAs with zeros for raster calculator
High_Rugosity_ras_calc <- High_Rugosity_ras
High_Rugosity_ras_calc[is.na(High_Rugosity_ras_calc[])] <- 0 

writeRaster(High_Rugosity_ras_calc, filename= "./Connectivity_between_MPA_Stephen/Rockfish_project/output_rasters/habitat_quality/inputs/High_rug.tiff", proj4string = High_Rugosity_ras@crs, overwrite=TRUE)

########################################################################
########################################################################
#Kelp beds
Kelp_beds <- readOGR("./Connectivity_between_MPA_Stephen/Rockfish_project/input_data/bcmca_eco_kelp_giantkelp_polygons_data", "BCMCA_ECO_Kelp_GiantKelp_Polygons_DATA")

#project
Kelp_beds <- spTransform(Kelp_beds, Ecozone_mask@proj4string)

#clip
dataframed <- proc.time() 
Kelp_beds <- gIntersection(Kelp_beds, Bathymetry, byid = FALSE) #couple seconds with my BCMCA layer
proc.time() - dataframed 

Kelp_beds_Value <- 1
Kelp_beds_Value <- as.data.frame(Kelp_beds_Value)

Kelp_beds <- SpatialPolygonsDataFrame(Kelp_beds, Kelp_beds_Value)
plot(Kelp_beds)

#Turn into raster
r <- raster(extent(Bathymetry))
projection(r) <- proj4string(Kelp_beds)
res(r) <- 1000 #in metres, size of each cell (e.g. 5000 = 5000x5000m) CHANGE BACK TO 5000

Kelp_beds_ras <- rasterize(Kelp_beds, r, field = "Kelp_beds_Value") #5 seconds

#Need to replace NAs with zeros for raster calculator
Kelp_beds_ras_calc <- Kelp_beds_ras
Kelp_beds_ras_calc[is.na(Kelp_beds_ras_calc[])] <- 0 

writeRaster(Kelp_beds_ras_calc, filename= "./Connectivity_between_MPA_Stephen/Rockfish_project/output_rasters/habitat_quality/inputs/Kelp_bed.tif", proj4string = Kelp_beds_ras@crs, overwrite=TRUE)



########################################################################
########################################################################
########################################################################
########################################################################
###[3] Raster calculator

Rockfish_habitat_suit <- raster::overlay(Benthic_ras_calc, High_Rugosity_ras_calc, Kelp_beds_ras_calc, fun = function(x,y,z){(x+y+z)})
plot(Rockfish_habitat_suit)

writeRaster(Rockfish_habitat_suit, filename= "./Connectivity_between_MPA_Stephen/Rockfish_project/output_rasters/habitat_quality/Rock_habsuit.tiff", proj4string = Kelp_beds_ras@crs, overwrite=TRUE)

proc.time() - full.run.time
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

