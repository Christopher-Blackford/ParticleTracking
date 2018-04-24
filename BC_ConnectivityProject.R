##############################################################################################################################
##############################################################################################################################
###BC_ConnectivityProject.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

###READ.ME:
#This file takes: 
#1) Input (release) and output (settle) data from the LTRANS model of larval dispersal
#2) A shapefile of the BC inshore area
#
#To build shapefiles showing connectivity on the BC coast
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


pld <- c(30)

year <- as.numeric(c(1998:2007)) 
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
### [2] Setting up study extent you will be using to clip your larval release points to your BC study extent

#Clipping to your study extent
BC_project_extent <- readOGR("./BC_ConnectivityProject/StudyExtent/Sarah_extent", "PU_ResizedTo10km")
My_BC_projection <- proj4string(BC_project_extent)
#Loading Remi's grid where larvae were released
grid <- readOGR("./cuke_present/StudyExtent/Starting_grid", "grid")
NAD_projection <- proj4string(grid)

#Dissolve into one polygon since so you can change grid dimensions
grid <- spTransform(grid, BC_project_extent@proj4string) #For some reason not "identical" to My_BC_projection, check later
grid <- gUnaryUnion(grid)

#Intersecting - don't know why this works and ConPoly2 <- grid[Ecozone_mask,] doesn't
ConPoly <- gIntersection(grid, BC_project_extent, byid = TRUE, drop_lower_td = TRUE) #This works, but you'll have to choose a shapefile that includes islands and doesn't cut-off at rivers 
#Adding dataframe so you can create a shapefile of new study extent
row.names(ConPoly) <- gsub(x = row.names(ConPoly), pattern = "1 ", ""); Poly_ID = row.names(ConPoly) 
temp <- as.data.frame(Poly_ID)
row.names(temp) <- row.names(ConPoly)
ConPoly <- SpatialPolygonsDataFrame(ConPoly, temp)

#Removing cells that weren't 95% of 10x10km
ConPoly@data$Area <- rgeos::gArea(ConPoly, byid = TRUE)/100000000
ConPoly <- ConPoly[(ConPoly$Area >= 0.95),]
plot(ConPoly)

writeOGR(ConPoly, dsn = "./BC_ConnectivityProject/StudyExtent/BCProject_Extent", layer = "BCProject_Extent", 
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

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
    #[4] Creating connectivity dataframes
    
    #Showing where each larvae begings and ends
    Release_df <- subset(Con_df, select = c(long0, lat0, Z0, larvae_ID, rday))
    Settle_df <- subset(Con_df, select = c(long, lat, Z, larvae_ID, year, rday))
    
    #Associate released points with where they were released from
    xy <- subset(Release_df, select = c(long0, lat0))
    
    Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
    Released_larvae <- spTransform(Released_larvae, ConPoly@proj4string) #use your custom BC projection for this
    #This shapefile ^^^ contains all the points in Remi's grid, including those outside of your study extent
    
    #Finding which polygons released larvae are in
    Released_larvae <- my.point.in.poly(Released_larvae, ConPoly) #takes many minutes
    #If you plot this file^^^, it only includes points within your study extent - it's doing some sort of merge with Conpoly
    writeOGR(Released_larvae, dsn = "./BC_ConnectivityProject/StudyExtent/Larvae_release", layer = "Larvae_release", driver = "ESRI Shapefile", 
             verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
    
    
    if (Make_depth_layers == TRUE & year_time == 1 & pld_time == 1){
      Released_dataframe <- Released_larvae@data
      Released_dataframe <- Released_dataframe[,c("Poly_ID", "Area")]
      counted <- count(Released_dataframe, Poly_ID)
      counted$larv <- counted$n/61
      Released_dataframe <- dplyr::distinct(Released_dataframe, Poly_ID, .keep_all = TRUE)
      
      Released_layer <- sp::merge(ConPoly, Released_dataframe, by = "Poly_ID", all.x = FALSE)
      
      writeOGR(Released_layer, dsn = "./BC_ConnectivityProject/StudyExtent/Larvae_release", layer = "Larvae_release_grid", driver = "ESRI Shapefile", 
               verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
      }
    
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
    Con_df <- merge(Released_larvae@data, Settled_larvae@data, by = "larvae_ID", all = T)
    #Need to remove larvae spawning outside of my study extent
    Con_df <- Con_df[complete.cases(Con_df[,"Poly_ID.x"]),]
    
    
    ########################################################################  
    ###[4b] Controlling for biased release
    #Need to set upper bound on how many larvae can be released from grid cell since areas closer to shore will have more larvae released
    #Don't need to set lower bound since I'm assuming cells where not a lot of larvae released are areas where not a lot of that depth class exists
    if (Bias_release_files_preloaded == FALSE & year_time == 1 & pld_time == 1){
      
      Biased_release <- Con_df
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
      write.csv(Biased_release, "./BC_ConnectivityProject/output/Release_bias/ReleaseBias.csv")
      
      ###when Bias_release_files_preloaded == TRUE
    } else{print(paste0("Loading in previous ", Habitat_classes_names[i], " larval release file")) 
      Biased_release <- read.csv("./BC_ConnectivityProject/output/Release_bias/ReleaseBias.csv")
    }

    ###End of biased release file creation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    #In cases where the above loop isn't true, make sure you have release dataframe to merge with habitat dataframe later
    ###Taking biased release file and un-biasing your depth class files
    Released_larvae_df <- merge(Biased_release, Con_df, by = c("larvae_ID", "long0", "lat0", "Z0", "Poly_ID.x")) #this is a dumb way to merge
    Released_larvae_df <- Released_larvae_df[(Released_larvae_df$Larv_code <= 100),] #change this 30 to min_release or something more automated at some point
    #Remove NAs for when settled and released don't line up
    Released_larvae_df <- Released_larvae_df[complete.cases(Released_larvae_df[,"Poly_ID.y"]),]
    Released_larvae_df <- Released_larvae_df[with(Released_larvae_df, order(Poly_ID.x, Poly_ID.y)), ]
    #write out final release bias file so it's obvious you only include 100 per cell
    write.csv(Released_larvae_df, paste0("./output_keep/release_settlement/hexagon/Release_bias/", Habitat_classes_names[i], "/merged_bias_file/", Habitat_classes_names[i], "_bias_controlled.csv"))
    
    
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
      
      
  #remi_in_poly <- my.point.in.poly(Released_larvae, grid)
    
    