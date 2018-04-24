#Trying to get boundary length modifier

#
##
###
####
#####
#Clear workspace
rm(list=ls())

require(rgdal)
require(rgeos)
require(dplyr)

my_resolution <- 10000 #defines raster cell size and controls for biased larvae release

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

#Adding dataframe so you can create a shapefile of new study extent
ConPoly_ID <- 1 
ConPoly <- SpatialPolygonsDataFrame(ConPoly, as.data.frame(ConPoly_ID))

#Turn your study extent into raster to analyse movement between cells
require(raster)
r <- raster(extent(ConPoly))
projection(r) <- proj4string(ConPoly)
res(r) <- my_resolution #in metres, size of each cell (e.g. 5000 = 5000x5000m)
ConGrid <- rasterize(ConPoly, r) #5 seconds

#Convert back to polygon for useable shapefile
ConPoly <- rasterToPolygons(ConGrid, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) #couple of seconds, can be much longer depending on file size

#Adding in unique IDs for each polygon
ConPoly@data$Poly_ID <- as.numeric(row.names(ConPoly))
#And removing useless "layer" column from attribute table
ConPoly <- ConPoly[,!(names(ConPoly) %in% "layer")]
rm(grid, ConPoly_ID, ConGrid, Ecozone_mask)
########################################################################
rm(r, My_BC_projection, NAD_projection)




#This should work
Boundary_length <- gTouches(ConPoly, byid = TRUE)
Boundary_length <- as.table(Boundary_length)
Boundary_length <- as.data.frame(Boundary_length)


Boundary_length <- Boundary_length[Boundary_length$Freq == "TRUE", ]
Boundary_length$Freq <- my_resolution
Boundary_length <- dplyr::rename(Boundary_length, id1 = Var1, id2 = Var2)
Boundary_length$Min_cell <- NA
Boundary_length$Max_cell <- NA

Boundary_length$id1 <- as.numeric(Boundary_length$id1)| Boundary_length$id2 <- as.numeric(Boundary_length$id2)
Boundary_length <- Boundary_length[with(Boundary_length, order(id1, id2)), ]
row.names(Boundary_length) <- 1:nrow(Boundary_length)

for (i in 1:nrow(Boundary_length)){
  Boundary_length$Min_cell[i] <- as.character(min(c(Boundary_length$id1[i], Boundary_length$id2[i])))
  Boundary_length$Max_cell[i] <- as.character(max(c(Boundary_length$id1[i], Boundary_length$id2[i])))
}

Boundary_length$Unique_cell <- paste0(Boundary_length$Min_cell, Boundary_length$Max_cell)
Boundary_length <- Boundary_length[!duplicated(Boundary_length$Unique_cell) ,]
Boundary_length <- Boundary_length[c("id1", "id2", "Freq")]





#Robert Bivand prefers this
require(spdep)
test2 <- spdep::poly2nb(ConPoly)
test2 <- nb2mat(test2, style = "B", zero.policy = TRUE)

test3 <- as.table(test2)

test3 <- as.data.frame(test2)
