###Getting total larvae release file
Out_release <- Con_df ; Out_release <- subset(Out_release, select = c(long0, lat0, Z0, larvae_ID))
xy <- subset(Out_release, select = c(long0, lat0))
Out_release <- SpatialPointsDataFrame(coords = xy, data = Out_release, proj4string = CRS(NAD_projection))
Out_release <- spTransform(Out_release, ConPoly@proj4string) #use your custom BC projection for this
writeOGR(Out_release, dsn = "./output_keep/release_settlement/zLarvae_release_locations/remi_location", layer = "BC_larvae_release", driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

