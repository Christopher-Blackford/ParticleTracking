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

#Clear workspace
rm(list=ls())

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
require(data.table); require(tidyverse); require(rgdal); require(rgeos)
require(maptools); require(sp); require(igraph);
#require(spatialEco);

###################Loading functions:
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/functions/my_point_in_poly.R")

###################Initialize run with these important parameters

Cell_cutoff_threshold <- 0 #Threshold determining percent the clipped cell needs to be compared to Sarah's extent to be included. 0 = All cells included. 1 = Only cells completely within Sarah's extent included.

Create_cutoff_graph <- FALSE #Do you want to create a histogram of how many cells get cut off by how much

Calculate_Larvae_per_cell <- FALSE #Do you want to create a csv showing how many larvae are released in each cell

Create_All_poly_pairs <- FALSE

pld <- c(30,60,120)

year <- as.numeric(c(1998:2007))
# ^ is equivalent to year <- c(1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007)

########################################################################
########################################################################
########################################################################
### [1] Loading up larval release points
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/1_Loading_up_larval_release_points.R")


########################################################################
########################################################################
########################################################################
### [2] Setting up study extent you will be using to clip your larval release points to your BC study extent

#Clipping to your study extent
BC_project_extent <- readOGR("./BC_ConnectivityProject/BC_StudyExtent/Sarah_extent", "PU_ResizedTo10km")
My_BC_projection <- proj4string(BC_project_extent)
#Loading Remi's grid where larvae were released
grid <- readOGR("./BC_ConnectivityProject/BC_StudyExtent/Remi_grid", "grid")
NAD_projection <- proj4string(grid)

#Dissolve into one polygon since so you can change grid dimensions
grid <- spTransform(grid, BC_project_extent@proj4string)
grid <- gUnaryUnion(grid)

#Intersecting Sarah and Remi's extent
ConPoly <- gIntersection(grid, BC_project_extent, byid = TRUE, drop_lower_td = TRUE) #This works, but you'll have to choose a shapefile that includes islands and doesn't cut-off at rivers 

#Adding dataframe so you can create a shapefile of new study extent
row.names(ConPoly) <- gsub(x = row.names(ConPoly), pattern = "1 ", ""); Poly_ID = row.names(ConPoly)
temp <- as.data.frame(Poly_ID); temp$Poly_ID <- as.numeric(as.character(temp$Poly_ID))+1; row.names(temp) <- row.names(ConPoly)
ConPoly <- SpatialPolygonsDataFrame(ConPoly, temp) #Poly_ID is equal to OBJECTID in Sarah's original file

#Finding area of each cell relative to Sarah's starting area
ConPoly@data$Area <- rgeos::gArea(ConPoly, byid = TRUE)/10^8 

#Create graph of number of cells cutoff by your clip if you want
if(Create_cutoff_graph == TRUE){source("./BC_ConnectivityProject/subcode/Cutoff graph.R")}

#Removing cells that weren't a certain percentage of 10x10km
ConPoly <- ConPoly[(ConPoly$Area >= Cell_cutoff_threshold),]
#writeOGR(ConPoly, dsn = "./BC_ConnectivityProject/BC_StudyExtent/BCProject_Extent", layer = paste0("BCProject_Extent_", Cell_cutoff_threshold), 
 #        driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

rm(BC_project_extent, grid, temp, Poly_ID) #If code doesn't work, this may be the culprit
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
    
    ###Reshaping dataset to take filename info and turning it into columns - This process takes a long time ~ 5 - 10 minutes
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
    #[4] Identifying locations of larval release and settlement
    
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
    #writeOGR(Released_larvae, dsn = "./BC_ConnectivityProject/BC_StudyExtent/Larvae_release", layer = "Larvae_release", driver = "ESRI Shapefile", 
     #        verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
    
    if (year_time == 1 & pld_time == 1 & Calculate_Larvae_per_cell == TRUE){
      Released_dataframe <- Released_larvae@data
      Released_dataframe <- Released_dataframe[,c("Poly_ID", "Area")]
      counted <- dplyr::count(Released_dataframe, Poly_ID)
      counted <- dplyr::rename(counted, Larvae_release_over_season = n)
      counted$Larvae_release_daily <- counted$Larvae_release_over_season/61 #Average should be 100 I think (400 for 20x20 = 100 for 10x10)
      write.csv(counted, "./BC_ConnectivityProject/BC_StudyExtent/Larvae_release/csv/Larvae_per_cell.csv", row.names = FALSE)
      
      png("./BC_ConnectivityProject/BC_StudyExtent/Larvae_release/csv/Larvae release per cell - daily.png", width = 1200, height = 800)
      hist(counted$Larvae_release_daily, breaks = seq(from = 0, to = 1200, by = 50), main = "Number of larvae released (daily) per cell", xlab = "Number of larvae released")
      abline(v = mean(counted$Larvae_release_daily), col = "red"); text(x = 800, y= 250, labels = paste0("Average = ", mean(counted$Larvae_release_daily)))
      dev.off()
      
      png("./BC_ConnectivityProject/BC_StudyExtent/Larvae_release/csv/Larvae release per cell - total.png", width = 1200, height = 800)
      hist(counted$Larvae_release_over_season, breaks = seq(from = 0, to = 70000, by = 7000), main = "Number of larvae released (total) per cell", xlab = "Number of larvae released")
      abline(v = mean(counted$Larvae_release_over_season), col = "red"); text(x = 50000, y= 500, labels = paste0("Average = ", mean(counted$Larvae_release_over_season)))
      dev.off()
      #Released_dataframe <- dplyr::distinct(Released_dataframe, Poly_ID, .keep_all = TRUE)
      #Released_layer <- sp::merge(ConPoly, Released_dataframe, by = "Poly_ID", all.x = FALSE)
      #writeOGR(Released_layer, dsn = "./BC_ConnectivityProject/BC_StudyExtent/Larvae_release", layer = "Larvae_release_grid", driver = "ESRI Shapefile", 
       #        verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
      }else{counted <- read.csv("./BC_ConnectivityProject/BC_StudyExtent/Larvae_release/csv/Larvae_per_cell.csv")}
    
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
    Released_larvae_df <- Con_df[complete.cases(Con_df[,"Poly_ID.y"]),]
    

    ########################################################################
    ########################################################################
    ########################################################################
    #[5] Creating connectivity metrics
    
      my_table <- table(Released_larvae_df$Poly_ID.x, Released_larvae_df$Poly_ID.y) #Counts number of larvae moving from Poly.X to Poly.y
      
      #As dataframe
      df <- as.data.frame(my_table)
      df$Var1 <- as.character(df$Var1); df$Var1 <- as.numeric(df$Var1); df$Var2 <- as.character(df$Var2); df$Var2 <- as.numeric(df$Var2)
      df <- dplyr::rename(df, Poly_ID_Release = Var1, Poly_ID_Settle = Var2)
      
      
      #Can probably remove these columns for now
      df <- base::merge(df, counted, by.x = "Poly_ID_Release", by.y = "Poly_ID") #Make sure number of larvae being released from counted is accurate to the Release_df
      df$Percent <- (df$Freq/(df$Larvae_release_over_season))*100

      #Writing out csv transforms from shapefile to the attribute table of that shapefile
      my_directory <- paste0("./BC_ConnectivityProject/BC_output/Con_df/pld", pld[pld_time])
      dir.create(my_directory)
      dir.create(paste0(my_directory, "/years"))
      write.csv(df, paste0(my_directory, "/years/Con_df_year", year[year_time], "pld", pld[pld_time], ".csv"), row.names = F)
      
  } #closing year loop
  
  print(paste0("Finished pld", pld[pld_time]))
} #closing pld loop

rm(my_table)

########################################################################
########################################################################
########################################################################
#[6] Merging connectivity dataframes across years to get: 1) Dataframes describing connectivity 2) Shapefiles describing connectivity 

#MIGHT HAVE TO INCORPORATE LONGER FILENAME TO INCLUDE CUTOFF THRESHOLD

for (pld_time in 1:length(pld)){
  
  filenames <- list.files(path= paste0("./BC_ConnectivityProject/BC_output/Con_df/pld", pld[pld_time], "/years"), pattern= ".csv", full.names=TRUE, recursive=T)
  
  # load all files into a list
  datalist <- lapply(filenames, read.csv)
  
  # set the names of the items in the list, so that you know which file it came from
  datalist <- setNames(datalist,filenames)
  
  # Merging the dataframe (rbind the list)
  dataset <- data.table::rbindlist(datalist, idcol="filename")
  dataset <- data.frame(dataset)
  dataset$filename <- NULL
  rm(datalist)
  
  #Averaging - think this is better than summing for now
  Con_df_All <- group_by(dataset, Poly_ID_Release, Poly_ID_Settle)
  Con_df_All <- dplyr::summarise(Con_df_All, mean(Freq), mean(Percent), mean(Larvae_release_over_season), mean(Larvae_release_daily), var(Freq), var(Percent))
  
  #Renaming columns to remove brackets and writing out dataframe
  colnames(Con_df_All) <- gsub("\\(", "_", colnames(Con_df_All)); colnames(Con_df_All) <- gsub("\\)", "", colnames(Con_df_All))
  write.csv(Con_df_All, paste0("./BC_ConnectivityProject/BC_output/Con_df/pld", pld[pld_time], "/Con_df_pld", pld[pld_time], ".csv"), row.names = FALSE)
  
  #Preparing ajacency matrix to calculate centrality scores
  #~~~~~~~
  Polys_total <- unique(c(Con_df_All$Poly_ID_Release, Con_df_All$Poly_ID_Settle))
  All_poly_pairs <- data.frame(col1 = double(), col2 = double())
  
  #Get all possible combinations between Poly release and Poly settle, then merge with existing table to make square adjacency matrix
  #~~~~~~~
  if (Create_All_poly_pairs == TRUE){
    for (j in 1:length(Polys_total)){
      for (i in 1:length(Polys_total)){
        temp_row <- c(Polys_total[j],Polys_total[i])
        All_poly_pairs <- rbind(temp_row, All_poly_pairs)
      }
      print(j/length(Polys_total))
    }
    rm(Polys_total, temp_row)
    All_poly_pairs$included <- 1
    colnames(All_poly_pairs)[1] <- "Poly_ID_Release"; colnames(All_poly_pairs)[2] <- "Poly_ID_Settle"
    write.csv(All_poly_pairs, paste0("./BC_ConnectivityProject/BC_output/All_poly_pairs/All poly dataframe_Cutoff", Cell_cutoff_threshold, ".csv"), row.names = FALSE)
    
  }else{All_poly_pairs <- read.csv(paste0("./BC_ConnectivityProject/BC_output/All_poly_pairs/All poly dataframe_Cutoff", Cell_cutoff_threshold ,".csv"))}
  #~~~~~~~
  
  ########################################################################
  ####Connectivity matrix for raw larval transport
  Con_df_table <- Con_df_All[c("Poly_ID_Release", "Poly_ID_Settle", "mean_Freq")]
  
  Con_df_table <- base::merge(Con_df_table, All_poly_pairs, by = c("Poly_ID_Release", "Poly_ID_Settle"), all.y=TRUE)
  Con_df_table$mean_Freq[is.na(Con_df_table$mean_Freq)] <- 0
  #Con_df_table$included <- NULL
  
  #Creating connectivity matrix
  Con_matrix <- t(matrix(Con_df_table$mean_Freq, nrow = length(unique(Con_df_table$Poly_ID_Release)), ncol = length(unique(Con_df_table$Poly_ID_Settle))))
  colnames(Con_matrix) <- unique(Con_df_table$Poly_ID_Release); rownames(Con_matrix) <- unique(Con_df_table$Poly_ID_Settle)
  #Con_dataframe <- as.data.frame(Con_matrix)
  
  EC_df <- eigen(t(Con_matrix), symmetric = FALSE)
  Contribution <- Re(EC_df$vectors[,1])
  
  num_neg <- 0
  for (i in 1:length(Contribution)){if (Contribution[i] < 0){num_neg <- num_neg+1}}#If your eigenvectors are all negative need to record this to make it positive
  
  if (num_neg == 0){#do nothing
  }else if(num_neg == length(Contribution)){Contribution <- abs(Contribution)
  }else{(print("Double check eigenvectors"))} #Making eigenvectors possitive if all negative. Only way to "stop" a process is to put this in a source file and designate that option
  
  my_graph <- igraph::graph_from_adjacency_matrix(Con_matrix, mode = list("directed"))#This is right! Directed means the matrix entry gives the strength of connection
  #central_score <- igraph::eigen_centrality(my_graph, directed = TRUE)#Here, need directed to take the value of the connection into account
  bet_score <- igraph::betweenness(my_graph, directed = TRUE)
  
  Con_df_table <- data.frame(Poly_ID = row.names(Con_matrix), Eigen_C = Contribution, Between_C = bet_score)
  row.names(Con_df_table) <- 1:nrow(Con_df_table)
  write.csv(Con_df_table, paste0("./BC_ConnectivityProject/BC_output/Con_df/pld", pld[pld_time], "/Con_Freq_pld", pld[pld_time], ".csv"), row.names = FALSE)
  
  #Writing out shapefile
  ConPoly_new <- sp::merge(ConPoly, Con_df_table, by = "Poly_ID")
  #ConPoly_new <- spatialEco::sp.na.omit(ConPoly_new, col.name = "mean_Local_retention", margin=1) #Can turn off to look at all cells that don't have any retention
  
  shapefile_directory <- paste0("./BC_ConnectivityProject/BC_output/shapefiles/pld", pld[pld_time])
  dir.create(shapefile_directory)
  writeOGR(ConPoly_new, dsn = shapefile_directory, layer = paste0("Con_Freq_pld", pld[pld_time]),
           driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
  
  
  ########################################################################
  ####Connectivity matrix for percent larval transport
  Con_df_table <- Con_df_All[c("Poly_ID_Release", "Poly_ID_Settle", "mean_Percent")]
  
  Con_df_table <- base::merge(Con_df_table, All_poly_pairs, by = c("Poly_ID_Release", "Poly_ID_Settle"), all.y=TRUE)
  Con_df_table$mean_Percent[is.na(Con_df_table$mean_Percent)] <- 0
  Con_df_table$mean_Percent <- Con_df_table$mean_Percent*100 #Need to do this so betweeness centrality will register. It's okay since what matters is the relative difference between cells
  
  #Creating connectivity matrix
  Con_matrix <- t(matrix(Con_df_table$mean_Percent, nrow = length(unique(Con_df_table$Poly_ID_Release)), ncol = length(unique(Con_df_table$Poly_ID_Settle))))
  colnames(Con_matrix) <- unique(Con_df_table$Poly_ID_Release); rownames(Con_matrix) <- unique(Con_df_table$Poly_ID_Settle)
  Con_dataframe <- as.data.frame(Con_matrix)
  
  EC_df <- eigen(t(Con_matrix), symmetric = FALSE)
  Contribution <- Re(EC_df$vectors[,1])
  
  num_neg <- 0
  for (i in 1:length(Contribution)){if (Contribution[i] < 0){num_neg <- num_neg+1}}#If your eigenvectors are all negative need to record this to make it positive
  
  if (num_neg == 0){#do nothing
  }else if(num_neg == length(Contribution)){Contribution <- abs(Contribution)
  }else{(print("Double check eigenvectors"))} #Making eigenvectors possitive if all negative. Only way to "stop" a process is to put this in a source file and designate that option
  
  my_graph <- igraph::graph_from_adjacency_matrix(Con_matrix, mode = list("directed"))#This is right! Directed means the matrix entry gives the strength of connection
  #central_score <- igraph::eigen_centrality(my_graph, directed = TRUE)#Here, need directed to take the value of the connection into account
  bet_score <- igraph::betweenness(my_graph, directed = TRUE)
  
  Con_df_table <- data.frame(Poly_ID = row.names(Con_matrix), Eigen_C = Contribution, Between_C = bet_score)
  row.names(Con_df_table) <- 1:nrow(Con_df_table)
  write.csv(Con_df_table, paste0("./BC_ConnectivityProject/BC_output/Con_df/pld", pld[pld_time], "/Con_Percent_pld", pld[pld_time], ".csv"), row.names = FALSE)
  
  #Writing out shapefile
  ConPoly_new <- sp::merge(ConPoly, Con_df_table, by = "Poly_ID")
  
  writeOGR(ConPoly_new, dsn = shapefile_directory, layer = paste0("Con_Percent_pld", pld[pld_time]),
           driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
  
}

writeLines("Finished everything! Yay, yay, yay, you the best! \nClick here for dog: \nhttps://media.giphy.com/media/l2JhO5yaMLa93hVeM/giphy.gif") 

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
