##############################################################################################################################
##############################################################################################################################
###Consistent depth class map

#
##
###
####
#####
#Clear workspace
rm(list=ls())
full.run.time <- proc.time() # 3.5 hours for pld 19 with Bias_release files already created

###################Loading required packages:
require(data.table)
require(tidyverse)
require(rgdal)
require(rgeos)
require(maptools)
require(igraph)
require(spatialEco)
require(sp)

###################Loading functions:
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/functions/my_point_in_poly.R")


###################Initialize run with these important parameters
my_resolution <- 10000 #defines raster cell size and controls for biased larvae release

#any values should be the same
pld <- c(22)
year <- c(1998)

pld_time <- 1 ; year_time <- 1

########################################################################
########################################################################
### [1] Loading up larval release points
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/1_Loading_up_larval_release_points.R")

########################################################################
########################################################################
#[2] Setting up study extent you will be using to clip your larval release points to your BC study extent
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_study_extent_you_will_be_using.R")

########################################################################
########################################################################
#[3a] Identifying settlement locations and linking to release locations

#####Initializing
memory.limit(size=15000) #need to manually increase memory limit from default to process across years/pld

#Acquiring files
filenames <- list.files(path=paste0("./cuke_present/ConData/G",year[year_time]), pattern=glob2rx(paste0("*para    1",formatC(pld[pld_time]+1, width = 3, format = "d", flag = "0"),"*")), full.names=TRUE,recursive=T)
# filenames <- list.files(path=paste0("F:/MPA_particles/output_keep/G",year), pattern=glob2rx(paste0("*para    1",formatC(pld+1, width = 3, format = "d", flag = "0"),"*")), full.names=TRUE,recursive=T) #REMI COMMENT

# load all files into a list, read_csv is much faster than read.csv
datalist <- lapply(filenames, read_csv,
                   col_names = c("long","lat","Z","Out","site"),
                   col_types = cols("d","d","d","i","i")
)

# set the names of the items in the list, so that you know which file it came from
datalist <- setNames(datalist,filenames)

# rbind the list
dataset <- rbindlist(datalist, idcol="filename")
dataset$site <- NA
rm(datalist)

###This process takes a long time ~ 5 - 10 minutes
#Reshaping dataset to take filename info and turning it into columns
dataset <- dataset %>%
  mutate(temp=substr(filename,24,nchar(filename))) %>%
  separate(temp,c("temp_type_year","rday","bin","time"),"/",convert=TRUE) %>% 
  separate(temp_type_year,c("type","year"),sep=1,convert=TRUE) %>% 
  mutate(time=as.integer(substr(time,9,13))-1001)

#Linking release locations to settlement locations based on bin
for(i in unique(dataset$bin)){
  x <- rl$bin==i
  y <- dataset$bin==i
  dataset$long0[y] <- rl$long0[x]
  dataset$lat0[y] <- rl$lat0[x]
  dataset$Z0[y] <- rl$Z0[x]
  dataset$delay[y] <- rl$delay[x]
  dataset$site0[y] <- rl$site0[x]
  print(paste(i,sum(x),sum(y),sum(is.na(dataset$long0)))) # this is just to show its working
}

#Add larvae IDs to dataset
Con_df <- dataset
Con_df <- subset(Con_df, select = c(long0, lat0, Z0, long, lat, Z, year, rday))
Con_df$larvae_ID <- row.names(Con_df)

#Now you can remove some large files but only if you want to!
rm(dataset,filenames,x,y,i)

########################################################################
########################################################################
########################################################################
########################################################################
#[4] Creating connectivity dataframes for each habitat class from Con_df dataframe
habitat.dataframe.time <- proc.time() #13 minutes

######Creating names for depth classes
Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")

######Creating depth classes based on range of depths
Intertidal_depth_limit <- c(15, -15)
Nearshore_depth_limit <- c(-15, -60)
Offshore_depth_limit <- c(-60, -250)

Depths_list <- c(Intertidal_depth_limit, Nearshore_depth_limit, Offshore_depth_limit)
Depths_list <- unique(Depths_list)

########################################################################
#####[4a] Linking larvae release to larvae settlement in dataframe

for (i in (1:length(Habitat_classes_names))){
  
  Habitat_classes <- Habitat_classes_names[i]
  
  #Release locations for habitat class
  Habitat_classes <- Con_df[(Con_df$Z0 < Depths_list[i] & Con_df$Z0 >= Depths_list[i+1]),]
  #Settlement locations for habitat class
  Habitat_classes <- Habitat_classes[(Habitat_classes$Z < Depths_list[i] & Habitat_classes$Z >= Depths_list[i+1]),]
  
  #Showing where each larvae begings and ends
  Release_df <- subset(Habitat_classes, select = c(long0, lat0, Z0, larvae_ID))
  Settle_df <- subset(Habitat_classes, select = c(long, lat, Z, larvae_ID, year, rday))
  
  #Associate released points with where they were released from
  xy <- subset(Release_df, select = c(long0, lat0))
  
  Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
  Released_larvae <- spTransform(Released_larvae, ConPoly@proj4string) #use your custom BC projection for this
  #This shapefile ^^^ contains all the points in Remi's grid, including those outside of your study extent
  
  #Write out release grids for Remi's extent
  #if (Bias_release_files_preloaded == FALSE & year_time == 1 & pld_time == 1){
  #writeOGR(Released_larvae, dsn = "./output_keep/release_settlement/zLarvae_release_locations/remi_location", layer = paste0("Remi_release_", Habitat_classes_names[i]), driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)}
  
  #Finding which polygons released larvae are in
  Released_larvae <- my.point.in.poly(Released_larvae, ConPoly) #takes many minutes
  
  #Making the depth layers from released points
  Depth_class_dataframe <- Released_larvae@data
  Depth_class_dataframe$Depth_class <- Habitat_classes_names[i]
  Depth_class_dataframe <- Depth_class_dataframe[,c("Poly_ID", "Depth_class")]
  Depth_class_dataframe <- dplyr::distinct(Depth_class_dataframe, Poly_ID, .keep_all = TRUE)
  
  Unif_depth_class <- sp::merge(ConPoly, Depth_class_dataframe, by = "Poly_ID", all.x = FALSE)
  
  writeOGR(Unif_depth_class, dsn = paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/", Habitat_classes_names[i]), 
           layer = paste0(Habitat_classes_names[i]), driver = "ESRI Shapefile", 
           verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
  
  #For next step
  assign(Habitat_classes_names[i], Unif_depth_class)
  }

