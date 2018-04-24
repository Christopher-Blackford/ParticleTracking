########################################################################
########################################################################
########################################################################
########################################################################
#[2] Setting up hexagon study extent you will be using to clip your larval release points to your BC study extent

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

writeOGR(ConPoly, dsn = "./cuke_present/StudyExtent/Inshore_extent", layer = "Ecozone_Remi_extent_merge", 
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)


#Creating the hexagon layer
ConPoly_hex <- ConPoly

hex_area <- my_resolution^2
dist_between_centroids <- sqrt(2*hex_area/sqrt(3))

ConPoly_hex = gBuffer(ConPoly_hex, width = 0)
HexPts <- spsample(ConPoly_hex, type="hexagonal", cellsize = dist_between_centroids, offset = c(0, 0))
HexPols <- HexPoints2SpatialPolygons(HexPts)
plot(HexPols)

row.names(HexPols) <- as.character(1:length(HexPols))
Poly_ID <- 1:length(HexPols)

HexPols <- SpatialPolygonsDataFrame(HexPols, as.data.frame(Poly_ID))
ConPoly <- HexPols

rm(grid, ConPoly_ID, Ecozone_mask, hex_area, dist_between_centroids, Poly_ID, ConPoly_hex, HexPts, HexPols)

writeOGR(ConPoly, dsn = "./cuke_present/StudyExtent/Inshore_extent_grid", layer = "Hexagon_grid", 
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

