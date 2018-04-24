#Creating map of BC habitats
rm(list=ls())

#Loading libraries
require(raster)
require(sf)
require(tidyverse)
require(marmap)


################################################################################################################################################
################################################################################################################################################
#Loading in shapefiles

#Datafiles
Coastal_zone <- st_read("cuke_present/StudyExtent/Inshore_extent/Inshore_extent.shp")
Coastal_zone_proj4string <- st_crs(Coastal_zone)$proj4string

BC_land <- st_read("K:/Christopher_PhD/CH1_MPA/Study_extent/Westcoast_shapefiles/From_statistics_Canada/Projected/British_Columbia_BCan.shp")

Intertidal <- st_read("./output_keep/shapefiles/Intertidal/pld20/Intertidal_pld20.shp")

Nearshore <- st_read("./output_keep/shapefiles/Nearshore/pld20/Nearshore_pld20.shp")

Offshore <- st_read("./output_keep/shapefiles/Offshore/pld20/Offshore_pld20.shp")

####################################
#Basemap

#Extent
lonmin <- -134.5; lonmax <- -122.17; latmin <- 47.75 ;latmax <- 55.5

#Showing benthic profile
BC_Ocean <- getNOAA.bathy(lonmin,lonmax,latmin,latmax,resolution=1,keep=TRUE)
plot.bathy(BC_Ocean,image=T,drawlabels=TRUE)

#plotting to view
blues <- colorRampPalette(c("darkblue", "dodgerblue1"))
greys <- colorRampPalette(c(grey(0.4),grey(0.99)))

#plot.bathy(BC_Ocean, image = TRUE, land = TRUE, n=0,
 #          bpal = list(c(0, max(BC_Ocean), greys(100)),
  #                     c(min(BC_Ocean), 0, blues(100))))

#Projecting data
BC_temp <- fortify.bathy(BC_Ocean)
BC_temp$z[BC_temp$z>0] <- NA

df <- data.frame(BC_temp)

#raster package
BC_ras <- rasterFromXYZ(df)  #Convert first two columns as lon-lat and third as value                
crs(BC_ras) <- "+proj=longlat +datum=NAD83"

BC_poly <- rasterToPolygons(BC_ras)

BC_background <- st_as_sf(BC_poly)

BC_background <- st_transform(BC_background, crs = Coastal_zone_proj4string)

#plot(BC_background["geometry"])
#plot(Coastal_zone["geometry"], add=T)
#plot(Intertidal["geometry"], add=T)

rm(BC_Ocean, BC_temp, BC_poly, BC_ras)

My_plot <- ggplot()+
  #geom_sf(data = BC_land, fill = "brown", colour = NA)+
  geom_sf(data = BC_background, aes(fill=z), colour = NA)+
  #geom_sf(data = Coastal_zone, fill = NA)+ don't really need the zone
  geom_sf(data = Offshore, fill = "gold", colour = NA)+
  geom_sf(data = Nearshore, fill = "forestgreen", colour = NA)+
  geom_sf(data = Intertidal, fill = "red2", colour = NA)+
  xlim(c(-235717.8, 504501.7))+
  ylim(c(905891, 1588793))+
  theme(panel.grid.major = element_line(colour = "transparent"),
      panel.background = element_rect("grey90")
      )

#My_plot

ggsave("K:/Christopher_PhD/temp/MainFigure.png", width = 10, height = 6)

write_sf(Intertidal, dsn = "K:/Christopher_PhD/temp", layer = "Intertidal_areas", driver = "ESRI Shapefile", delete_layer = TRUE, update = TRUE)
write_sf(Nearshore, dsn = "K:/Christopher_PhD/temp", layer = "Nearshore_areas", driver = "ESRI Shapefile", delete_layer = TRUE, update = TRUE)
write_sf(Offshore, dsn = "K:/Christopher_PhD/temp", layer = "Offshore_areas", driver = "ESRI Shapefile", delete_layer = TRUE, update = TRUE)


#Backup xlim and ylim 
#back1 <- c(-235717.8, 504501.7)
#back2 <- c(905891, 1638793)


#ggplot()+
 # geom_sf(data = BC_land, fill = "brown", colour = NA)+
  #geom_sf(data = Coastal_zone, fill = NA)+
  #xlim(c(-235717.8, 504501.7))+
  #ylim(c(905891, 1638793))
