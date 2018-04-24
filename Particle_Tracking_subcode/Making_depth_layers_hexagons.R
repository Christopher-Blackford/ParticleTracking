#Making_depth_layers

#Making the depth layers from released points
Depth_class_dataframe <- Released_larvae@data
Depth_class_dataframe$Depth_class <- Habitat_classes_names[i]
Depth_class_dataframe <- Depth_class_dataframe[,c("Poly_ID", "Depth_class")]
Depth_class_dataframe <- dplyr::distinct(Depth_class_dataframe, Poly_ID, .keep_all = TRUE)

Unif_depth_class <- sp::merge(ConPoly, Depth_class_dataframe, by = "Poly_ID", all.x = FALSE)

writeOGR(Unif_depth_class, dsn = paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/", Habitat_classes_names[i]), 
         layer = paste0(Habitat_classes_names[i]), driver = "ESRI Shapefile", 
         verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

