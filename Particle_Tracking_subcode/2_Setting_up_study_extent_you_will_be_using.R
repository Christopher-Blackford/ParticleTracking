########################################################################
########################################################################
########################################################################
########################################################################
#[2] Setting up study extent you will be using to clip your larval release points to your BC study extent

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
head(ConPoly@data)
rm(grid, ConPoly_ID, ConGrid, Ecozone_mask)

#To get you a shapefile of your study extent for chapter 1
#For ConPoly with grids
writeOGR(ConPoly, dsn = "./output_keep/shapefiles", layer = paste0("Release_grid_", my_resolution/1000, "km"),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

#Without grid
#CH1_studyextent <- gUnaryUnion(ConPoly)
#StudyExtent <- 1
#CH1_studyextent <- SpatialPolygonsDataFrame(CH1_studyextent, as.data.frame(StudyExtent))
#writeOGR(CH1_studyextent, dsn = "K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Coastal_region", layer = "CH1_StudyExtent",
#        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
