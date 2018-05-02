##############################################################################################################################
##############################################################################################################################
###Present_particle_locations.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

###READ.ME:
#This file takes: 
#1) Input (release) and output (settle) data from the LTRANS model of larval dispersal
#2) A shapefile of the BC inshore area
#
#To build shapefiles showing  connectivity in the BC inshore region between areas
#with the same depth class.
#The analysis can be run across multiple years and for multiple PLD values

#
##
###
####
#####
#Clear workspace
rm(list=ls())
full.run.time <- proc.time() # 3.5 hours for pld 19 with Bias_release files already created

###################TABLE OF CONTENTS
###[1] Loading up larval release points
###[2] Setting up study extent you will be using to clip your larval release points to your BC study extent
###[3] Identifying settlement locations and linking to release locations
#Begin big loop
###[4] Creating connectivity dataframes for each habitat class from Con_df dataframe
#####[4a] Linking larvae release to larvae settlement in dataframe
###[#[4b] Controlling for biased release
###[5] Creating connectivity metrics for each depth class
#End big loop
###[6] Merging connectivity dataframes across years to get: 1) Dataframes describing connectivity 2) Shapefiles describing connectivity 
###

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

Bias_release_files_preloaded <- TRUE #TRUE = Will use prior data to control for bias in number of larvae release per cell. 
                                     #FALSE = Performs operation that randomly removes larval from polygons where too many larvae are released.

Make_count_csv <- FALSE
Make_depth_layers <- FALSE

my_resolution <- 10000 #defines hexagon cell size

pld <- c(22,56,48)
#Intertidal pld average: 22 #Nearshore pld average: 56 #Offshore pld average: 48
pld <- c(3, 4, 6, 9, 21, 24, 27, 28, 29, 31, 38, 40, 45, 52, 58, 60, 61, 78, 90, 91, 95, 105, 109, 120) #PLD for all my species

#pld 61 for year 2003 throws an error
pld <- c(61)
year <- as.numeric(c(2003)) 
# ^ is equivalent to year <- c(1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007)


########################################################################
########################################################################
########################################################################
########################################################################
### [1] Loading up larval release points

source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/1_Loading_up_larval_release_points.R")

########################################################################
########################################################################
########################################################################
########################################################################
#[2] Setting up study extent you will be using to clip your larval release points to your BC study extent

source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")

########################################################################
########################################################################
########################################################################
########################################################################
#[3a] Identifying settlement locations and linking to release locations
pld_time <- 1
year_time <- 1  

#####Initializing
memory.limit(size=15000) #need to manually increase memory limit from default to process across years/pld

for (pld_time in 1:length(pld)){
  
  for (year_time in 1:length(year)){
        
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
    
    #Gettting total larvae release
    #source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/Getting total larvae release file.R")
    
    #Now you can remove some large files but only if you want to!
    rm(dataset,filenames,x,y,i)

    ########################################################################
    ########################################################################
    ########################################################################
    ########################################################################
    #[4] Creating connectivity dataframes for each habitat class from Con_df dataframe

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
      #If you plot this file^^^, it only includes points within your study extent - it's doing some sort of merge with Conpoly
      
      if (Make_depth_layers == TRUE & year_time == 1 & pld_time == 1){
        source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/Making_depth_layers_hexagons.R")}

      if (Make_count_csv == TRUE & year_time == 1 & pld_time == 1){
        Released_dataframe <- Released_larvae@data
        counted <- count(Released_dataframe, Poly_ID)
        counted$larv <- counted$n/61
        write.csv(counted, paste0("./output_keep/release_settlement/zLarvae_release_locations/", Habitat_classes_names[i], "_counted.csv"))
        assign(paste0(Habitat_classes_names[i], "_counted"), counted)}
      
      #Write out release grids for my (BC) extent
      #if (Bias_release_files_preloaded == FALSE & year_time == 1 & pld_time == 1){
        #writeOGR(Released_larvae, dsn = "./output_keep/release_settlement/zLarvae_release_locations/my_study_location", layer = paste0("BC_release_", Habitat_classes_names[i]), driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)}
  
      #Associate settled points with where they settled 
      xy <- subset(Settle_df, select = c(long, lat))
      
      Settled_larvae <- SpatialPointsDataFrame(coords = xy, data = Settle_df, proj4string = CRS(NAD_projection))
      Settled_larvae <- spTransform(Settled_larvae, ConPoly@proj4string) #use your custom BC projection for this
      
      #Finding which polygons settled larvae are in
      Settled_larvae <- my.point.in.poly(Settled_larvae, ConPoly) #takes many minutes
      
      #Join dataframes to make precursor to connectivity matrices
      Habitat_classes <- merge(Released_larvae@data, Settled_larvae@data, by = "larvae_ID", all = T)
      #Need to remove larvae spawning outside of my study extent
      Habitat_classes <- Habitat_classes[complete.cases(Habitat_classes[,"Poly_ID.x"]),]
      
      #For next step
      assign(Habitat_classes_names[i], Habitat_classes)
    }
    ########################################################################
    ########################################################################
    #Write out intertidal, nearshore, and offshore extents
    
    
    ########################################################################  
    ###[4b] Controlling for biased release
    #Need to set upper bound on how many larvae can be released from grid cell since areas closer to shore will have more larvae released
    #Don't need to set lower bound since I'm assuming cells where not a lot of larvae released are areas where not a lot of that depth class exists
    
    Habitat_classes <- list(Intertidal, Nearshore, Offshore)
    
    ###Creating biased release file if needed~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:length(Habitat_classes)){
    if (Bias_release_files_preloaded == FALSE & year_time == 1 & pld_time == 1){
        
    Biased_release <- Habitat_classes[[i]]
    set.seed(1)
    Biased_release$RAND <- sample(1:nrow(Biased_release), nrow(Biased_release), replace=F)
    Biased_release <- Biased_release[with(Biased_release, order(Poly_ID.x, RAND)), ]
      
    Biased_release$Larv_code <- NA #Need to fill it with something
    Biased_release$Larv_code[1] <- 1
    row.names(Biased_release) <- 1:nrow(Biased_release)
    
    ###Counting up number of larvae released in each cell. Takes over a day to do this for all habitat classes
    for (release_count in 2:nrow(Biased_release)){
      if (Biased_release$Poly_ID.x[release_count] == Biased_release$Poly_ID.x[release_count-1]){
        Biased_release$Larv_code[release_count] = Biased_release$Larv_code[release_count-1] + 1}
      else {Biased_release$Larv_code[release_count] = 1}
        }
  
    Biased_release <- Biased_release[c("larvae_ID", "long0", "lat0", "Z0", "Poly_ID.x", "RAND", "Larv_code")]
    write.csv(Biased_release, paste0("./output_keep/release_settlement/hexagon/Release_bias/", Habitat_classes_names[i], "/", Habitat_classes_names[i], "_61Release.csv"))
    assign(paste0(Habitat_classes_names[i], "_Release"), Biased_release)
    
    ###when Bias_release_files_preloaded == TRUE
    } else{print(paste0("Loading in previous ", Habitat_classes_names[i], " larval release file")) 
    Biased_release <- read.csv(paste0("./output_keep/release_settlement/hexagon/Release_bias/", Habitat_classes_names[i], "/", Habitat_classes_names[i], "_61Release.csv"))
    assign(paste0(Habitat_classes_names[i], "_Release"), Biased_release)
    }
      
    } #getting for Intertidal, Nearshore, and Offshore
    
    ###End of biased release file creation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    #In cases where the above loop isn't true, make sure you have release dataframe to merge with habitat dataframe later
    ###Taking biased release file and un-biasing your depth class files
    Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")
    Habitat_classes <- list(Intertidal, Nearshore, Offshore)
    Release_classes <- list(Intertidal_Release, Nearshore_Release, Offshore_Release)
    
    for (i in 1:length(Release_classes)){
      
      Released_larvae_df <- Release_classes[[i]]
      Released_larvae_df <- merge(Released_larvae_df, Habitat_classes[[i]], by = c("larvae_ID", "long0", "lat0", "Z0", "Poly_ID.x")) #this is a dumb way to merge
      Released_larvae_df <- Released_larvae_df[(Released_larvae_df$Larv_code <= 6100),] #old way was 100 #100 because 10x10km will have 100 larvae spawn. 61 because we are releasing these 100 larvae for 61 days
      #Remove NAs for when settled and released don't line up
      Released_larvae_df <- Released_larvae_df[complete.cases(Released_larvae_df[,"Poly_ID.y"]),]
      Released_larvae_df <- Released_larvae_df[with(Released_larvae_df, order(Poly_ID.x, Poly_ID.y)), ]
      #write out final release bias file so it's obvious you only include 100 per cell
      write.csv(Released_larvae_df, paste0("./output_keep/release_settlement/hexagon/Release_bias/", Habitat_classes_names[i], "/merged_bias_file/", Habitat_classes_names[i], "_61bias_controlled.csv"))
      
      assign(Habitat_classes_names[i], Released_larvae_df)
    }
    
    rm(Intertidal_depth_limit, Nearshore_depth_limit, Offshore_depth_limit, Habitat_classes_names, Habitat_classes, Depths_list)
    rm(xy, Release_df, Settle_df, Settled_larvae)
    rm(Biased_release, Intertidal_Release, Nearshore_Release, Offshore_Release) #Bias files
    ###now you have controlled for biased release~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ########################################################################
    ########################################################################
    
    
    ########################################################################
    ########################################################################
    ########################################################################
    ########################################################################
    #[5] Creating connectivity metrics for each depth class
    
    ########################################################################
    ###For looping across depth classes
    Habitat_classes <- list(Intertidal, Nearshore, Offshore)
    Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")
    
    for (i in 1:length(Habitat_classes_names)){
      
      my_table <- table(Habitat_classes[[i]]$Poly_ID.x, Habitat_classes[[i]]$Poly_ID.y)
      
      #As dataframe
      df <- as.data.frame(my_table)
      df$Var1 <- as.character(df$Var1)
      df$Var1 <- as.numeric(df$Var1)
      df$Var2 <- as.character(df$Var2)
      df$Var2 <- as.numeric(df$Var2)
      df <- dplyr::rename(df, Poly_ID_Release = Var1, Poly_ID_Settle = Var2)
      
      ###Local retention
      df_Local_Retention <- df[which(df$Poly_ID_Release == df$Poly_ID_Settle & df$Freq > 0),]
      df_Local_Retention <- df_Local_Retention[c("Poly_ID_Release", "Freq")]
      df_Local_Retention <- dplyr::rename(df_Local_Retention, Poly_ID = Poly_ID_Release, Local_retention = Freq)
      
      ###Eigenvector centrality
      df_Eigen_Centrality <- graph.data.frame(df)
      df_Eigen_Centrality <- eigen_centrality(df_Eigen_Centrality, directed=TRUE, weights=E(df_Eigen_Centrality)$Freq)$vector
      df_Eigen_Centrality <- as.data.frame(df_Eigen_Centrality)
      df_Eigen_Centrality$Poly_ID <- row.names(df_Eigen_Centrality)
      row.names(df_Eigen_Centrality) <- 1:nrow(df_Eigen_Centrality)
      df_Eigen_Centrality <- dplyr::rename(df_Eigen_Centrality, Eigenvector_centrality = df_Eigen_Centrality)
      
      #Merging these dataframes
      Released_per_gridcell <- Habitat_classes[[i]]
      Released_per_gridcell <- dplyr::count(Released_per_gridcell, Poly_ID.x)
      Released_per_gridcell <- dplyr::rename(Released_per_gridcell, Poly_ID = Poly_ID.x, Larvae_released = n)
      df_out <- merge(Released_per_gridcell, df_Local_Retention, by = "Poly_ID")
      df_out <- merge(df_out, df_Eigen_Centrality, by = "Poly_ID")
      df_out$percent_LR <- df_out$Local_retention/df_out$Larvae_released
      
      assign(Habitat_classes_names[i], df_out)
    }
    
    #Writing out csv transforms from shapefile to the attribute table of that shapefile
    
    Habitat_classes <- list(Intertidal, Nearshore, Offshore)
    Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")
    
    for (i in 1:length(Habitat_classes_names)){
      my_directory <- paste0("./output_keep/Con_df/hexagon/", Habitat_classes_names[i], "/pld", pld[pld_time])
      dir.create(my_directory)
      dir.create(paste0(my_directory, "/years"))
      write.csv(Habitat_classes[[i]], paste0(my_directory, "/years/", Habitat_classes_names[i], "_year", year[year_time], "pld", pld[pld_time], ".csv"), row.names = F)
    }
  
  } #closing year loop
  
  print(paste0("Finished pld", pld[pld_time]))
} #closing pld loop

rm(my_table)

########################################################################
########################################################################
########################################################################
########################################################################
#[6] Merging connectivity dataframes across years to get: 1) Dataframes describing connectivity 2) Shapefiles describing connectivity 

Habitat_classes_names <- c("Intertidal", "Nearshore", "Offshore")

for (pld_time in 1:length(pld)){
  
  for (j in 1:length(Habitat_classes_names)){
    
    filenames <- list.files(path= paste0("./output_keep/Con_df/hexagon/", Habitat_classes_names[j], "/pld", pld[pld_time], "/years"), pattern= ".csv", full.names=TRUE, recursive=T)
    
    # load all files into a list
    datalist <- lapply(filenames, read.csv)
    
    # set the names of the items in the list, so that you know which file it came from
    datalist <- setNames(datalist,filenames)
    
    # Merging the dataframe (rbind the list)
    dataset <- data.table::rbindlist(datalist, idcol="filename")
    dataset <- data.frame(dataset)
    dataset$filename <- NULL
    rm(datalist)
    
    #Averaging
    Con_df <- group_by(dataset, Poly_ID)
    Con_df <- dplyr::summarise(Con_df, mean(Larvae_released), mean(Local_retention), sd(Local_retention), mean(Eigenvector_centrality), sd(Eigenvector_centrality), mean(percent_LR))
    
    #Renaming and writing out csv
    colnames(Con_df) <- gsub("\\(", "_", colnames(Con_df)); colnames(Con_df) <- gsub("\\)", "", colnames(Con_df))
    write.csv(Con_df, paste0("./output_keep/Con_df/hexagon/", Habitat_classes_names[j], "/pld", pld[pld_time], "/", Habitat_classes_names[j], "_pld", pld[pld_time], ".csv"), row.names = FALSE)
    
    #Writing out shapefile
    ConPoly_new <- sp::merge(ConPoly, Con_df, by = "Poly_ID")
    ConPoly_new <- spatialEco::sp.na.omit(ConPoly_new, col.name = "mean_Local_retention", margin=1) #Can turn off to look at all cells that don't have any retention
    
    shapefile_directory <- paste0("./output_keep/shapefiles/hexagon/", Habitat_classes_names[j], "/pld", pld[pld_time])
    dir.create(shapefile_directory)
    writeOGR(ConPoly_new, dsn = shapefile_directory, layer = paste0(Habitat_classes_names[j], "_pld", pld[pld_time]),
             driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
  }
  
}

rm(i, j)

writeLines("Finished everything! Yay, yay, yay, you the best! \nClick here for dog: \nhttps://media.giphy.com/media/l2JhO5yaMLa93hVeM/giphy.gif") 
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


########################################################################
########################################################################


#Write out intertidal, nearshore, and offshore spatialextents
#Poly_ID_release <- unique(na.omit(Habitat_classes$Poly_ID.x))
#Poly_ID_settle <- unique(na.omit(Habitat_classes$Poly_ID.y))
#Study_extent <- as.data.frame(unique(append(Poly_ID_release, Poly_ID_settle))); rm(Poly_ID_release, Poly_ID_settle)
#names(Study_extent)[1] <- "Poly_ID"
#Study_extent <- dplyr::arrange(Study_extent, Poly_ID) #fixing rownames
#Study_extent <- merge(ConPoly, Study_extent, by = "Poly_ID", all.x = FALSE)
#writeOGR(Study_extent, dsn = paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/", Habitat_classes_names[i]), layer = paste0(Habitat_classes_names[i], "_", my_resolution/1000, "km"), 
 #        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
