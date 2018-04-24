#############################################################################################################################
##############################################################################################################################
###Connectivity_between_MPA.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)

###READ.ME:
#This file takes: 
#1) Input (release) and output (settle) data from the LTRANS model of larval dispersal
#2) A shapefile of the BC inshore area
#3) A shapefile of current MPA locations in BC
#
#To build shapefiles showing connectivity in the BC inshore region between 
#current distribution of MPAs
#The analysis can be run across multiple years and for multiple PLD values

#
##
###
####
#####

#Clear workspace
rm(list=ls())
full.run.time <- proc.time() # 33 minutes

###################TABLE OF CONTENTS
###[1] Loading up larval release points
###[2]Choosing year of release and pld you are tracking
###[3] Identifying settlement locations and linking to release locations
###[4] Setting up study extent you will be using to clip your larval release points to your BC study extent
###[5] Removing some of the MPAs that got sliced to thin clipping to Remi's extent
###[6] Creating dataframe describing release and settlement of each particle
###[7] Creating connectivity tables - (down) column donates to (across) row
###[8] Creating shapefiles with larval dispersal data represented

###################Loading required packages:
require(plyr)
require(data.table)
require(tidyverse)
require(readr)
require(rgdal)
require(rgeos)
require(maptools)
require(spatialEco)

###################Custom functions:
#Remove NAs in sp Dataframe object
#   x           sp spatial DataFrame object
#   margin      Remove rows (1) or columns (2) 
sp.na.omit <- function(x, margin=1) {
  if (!inherits(x, "SpatialPointsDataFrame") & !inherits(x, "SpatialPolygonsDataFrame")) 
    stop("MUST BE sp SpatialPointsDataFrame OR SpatialPolygonsDataFrame CLASS OBJECT") 
  na.index <- unique(as.data.frame(which(is.na(x@data),arr.ind=TRUE))[,margin])
  if(margin == 1) {  
    cat("DELETING ROWS: ", na.index, "\n") 
    return( x[-na.index,]  ) 
  }
  if(margin == 2) {  
    cat("DELETING COLUMNS: ", na.index, "\n") 
    return( x[,-na.index]  ) 
  }
}

########################################################################
########################################################################
########################################################################
########################################################################
### [i] Looping across multiple size_reduction_thresholds

size_reduction_threshold <- c(85, 90, 95, 999)

length(size_reduction_threshold)

for (size_reduction_loop in 1:length(size_reduction_threshold)){
  
########################################################################
########################################################################
########################################################################
########################################################################
### [1] Loading up larval release points

#Acquiring files
filenames <- list.files(path = "./cuke_present/ReleaseLocations", pattern="rl_.", full.names=TRUE,recursive=T)

# load all files into a list, read_csv is much faster than read.csv
rllist <- lapply(filenames, read_csv,
                 col_names = c("long0","lat0","Z0","delay","site0"),
                 col_types = cols("d","d","i","i","i")
)

# set the names of the items in the list, so that you know which file it came from
rllist <- setNames(rllist,filenames)

# rbind the list
rl <- rbindlist(rllist, idcol="filename")

rl$bin <- as.numeric(gsub(".*rl_|.txt.*", "",rl$filename))
head(rl)
rm(rllist, filenames)

#Creating csv file ith all starting locations
#write.csv(rl, file="./output/release_settlement/Remi_release_lat_long.csv", row.names = F)


########################################################################
########################################################################
########################################################################
########################################################################
###[2] Setting up study extent you will be using to clip your larval release points to your BC study extent

#All MPAs in the region 
All_MPAs <- readOGR("./Connectivity_between_MPA_Sarah/study_extent/raw_MPA", "NSBMPAsForChris_16Jun17")
#Correct for updated shapefile with different MPA FID names
All_MPAs@data$FID_NSB_MP <- All_MPAs@data$NSB_MPA_ID
row.names(All_MPAs@data) <- All_MPAs@data$FID_NSB_MP #Change row.names because FID starts at 0 and you want it to start at 1
head(All_MPAs@data)
BC_projection <- proj4string(All_MPAs)

#Loading Remi's grid where larvae were released
grid <- readOGR("./cuke_present/StudyExtent/Starting_grid", "grid")
NAD_projection <- proj4string(grid)
proj4string(grid)
#Dissolve into one polygon since so you can change grid dimensions
grid <- spTransform(grid, BC_projection)
grid <- gUnaryUnion(grid)

#Loading MaPP boundary layer
MaPP_extent <- readOGR("./Connectivity_between_MPA_Sarah/study_extent/MaPP_Marine_Area_EstuaryCorrected", "MaPP_Marine_Area_EstuaryCorrected")
MaPP_extent@data

#Intersecting - don't know why this works and ConPoly2 <- grid[Ecozone_mask,] doesn't
MPA_mask <- gBuffer(All_MPAs, byid=TRUE, width=0) #Need to do this to avoid ring intersection
row.names(MPA_mask@data) <- MPA_mask@data$FID_NSB_MP #Change row.names because FID starts at 0 and you want it to start at 1
head(MPA_mask@data)
MPA_mask_id <- as.character(MPA_mask@data$FID_NSB_MP)
ConPoly <- gIntersection(grid, MPA_mask, byid = TRUE, id = MPA_mask_id) 
#Clip to MaPP extent (unnecessary?)
ConPoly <- gIntersection(MaPP_extent, ConPoly, byid = TRUE, id = row.names(ConPoly)) 


#Adding dataframe so you can create a shapefile of new study extent
#Clipped dataframe
ConPoly_ID <- row.names(ConPoly)
ConPoly_ID <- as.numeric(ConPoly_ID)
ConPoly_ID <- as.data.frame(ConPoly_ID)
row.names(ConPoly_ID) <- ConPoly_ID$ConPoly_ID


#Original dataframe
#ConPoly_data <- as.data.frame(MPA_mask[ConPoly_ID$ConPoly_ID, ]) Not sure why this was used earlier... delete in a week or so
ConPoly_data <- MPA_mask@data

MPAS <- SpatialPolygonsDataFrame(ConPoly, ConPoly_ID)
MPAS@data <- plyr::rename(MPAS@data, c("ConPoly_ID" = "FID_NSB_MP"))

#Write out layer
#writeOGR(MPAS, dsn = "./Connectivity_between_MPA_Sarah/output_keep/shapefiles", layer = "MPA_clipped",
         #driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

rm(grid, ConPoly_ID, ConPoly, ConPoly_data, MPA_mask_id, MPA_mask)

########################################################################
###Removing some of the MPAs that got sliced to thin clipping to Remi's extent

##size_reduction_threshold <- 10 #You are only left with 9 if you go down to 0

#Compare how much smaller clipped MPA layer is to All_MPAs file
All_MPAs@data$Merged_area <- gArea(All_MPAs, byid = TRUE)

MPA_clipped_size <- MPAS
MPA_clipped_size@data$Clip_Area <- gArea(MPA_clipped_size, byid = TRUE)

MPA_clipped_size <- sp::merge(MPA_clipped_size, All_MPAs@data, by = "FID_NSB_MP")
MPA_clipped_size@data$size_reduction <- 100*(1 - MPA_clipped_size@data$Clip_Area/MPA_clipped_size@data$Merged_area)
row.names(MPA_clipped_size@data) <- MPA_clipped_size@data$FID_NSB_MP #Change row.names 

Size_reduction_df <- MPA_clipped_size@data

#write.csv(Size_reduction_df, "./Connectivity_between_MPA_Sarah/output_keep/size_reduction.csv", row.names = F)

#Histogram of how many MPAs got clipped by how much
Size_histogram <- ggplot(Size_reduction_df, aes(size_reduction)) +
  geom_histogram(binwidth = 5, fill = "#FFFFFF", colour = "black") +
  labs(title = "Histogram of MPA loss", x = "Percent loss", y = "Count") +
  theme(
    plot.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.line = element_line("black"),
    panel.background = element_blank()
  )
Size_histogram

#Removing MPAs that were too clipped by Remi's extent (based on percent loss)
MPAS <- MPA_clipped_size[MPA_clipped_size@data$size_reduction <= size_reduction_threshold[size_reduction_loop],]
MPAS_loop <- MPAS

########################################################################
########################################################################
########################################################################
########################################################################
### [3] Choosing year of release and pld you are tracking

memory.limit(size=15000)

# List the particle tracking files for that particular year and pld
year <- c(1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007)
pld <- c(30, 60, 120)

###To Loop: Remove hastags below and add closing brackets, delete temporary year_time and pld_time, and concatate year and pld vectors above
###
for (year_time in 1:length(year)){
  
  for (pld_time in 1:length(pld)){

###temporary year_time and pld_time for non loops
#year_time <- 1
#pld_time <- 1
#i <- 1
########################################################################
########################################################################
########################################################################
###[4] Identifying settlement locations and linking to release locations

#Acquiring files
filenames <- list.files(path=paste0("./cuke_present/ConData/G", year[year_time]), pattern=glob2rx(paste0("*para    1",formatC(pld[pld_time]+1, width = 3, format = "d", flag = "0"),"*")), full.names=TRUE,recursive=T)

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
  # mutate(temp=substr(filename,25,nchar(filename))) %>% # you probably want this back to 24? #REMI COMMENT
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
rm(filenames,x,y,i)

#Write out dataset - takes a long time
#write.csv(dataset, "./output/connectivity_tables/dataset.csv", row.names = F) #6 minutes

#Add larvae IDs to dataset
Con_df <- dataset
Con_df <- subset(Con_df, select = c(long0, lat0, Z0, long, lat, Z, year, rday))
Con_df$larvae_ID <- row.names(Con_df)

#Write out Con_df - takes a long time
#write.csv(Con_df, "./output/connectivity_tables/Con_df.csv", row.names = F)

#Now you can remove some large files but only if you want to!
rm(dataset)


########################################################################
########################################################################
########################################################################
########################################################################
###[5] Creating dataframe describing release and settlement of each particle
#Clipping to FID_NSB_MP to do points in poly
MPAS <- MPAS_loop[,"FID_NSB_MP"]

MPA.dataframe.time <- proc.time() #6 minutes

#####Showing where each larvae begings and ends
Release_df <- subset(Con_df, select = c(long0, lat0, Z0, larvae_ID))
Settle_df <- subset(Con_df, select = c(long, lat, Z, larvae_ID, year, rday))
rm(Con_df) #to free up space

#Associate released points with where they were released from
xy <- subset(Release_df, select = c(long0, lat0))
Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
Released_larvae <- spTransform(Released_larvae, MPAS@proj4string) #use your custom BC projection for this
Released_larvae <- Released_larvae[MPAS,]
#Finding which polygons released larvae are in
Released_larvae <- point.in.poly(Released_larvae, MPAS) #takes many minutes


#Associate settled points with where they settled 
xy <- subset(Settle_df, select = c(long, lat))
Settled_larvae <- SpatialPointsDataFrame(coords = xy, data = Settle_df, proj4string = CRS(NAD_projection))
Settled_larvae <- spTransform(Settled_larvae, MPAS@proj4string) #use your custom BC projection for this
Settled_larvae <- Settled_larvae[MPAS,]
#Finding which polygons settled larvae are in
Settled_larvae <- point.in.poly(Settled_larvae, MPAS) #takes many minutes


#Join dataframes to make precursor to connectivity matrices
MPA_df <- merge(Released_larvae@data, Settled_larvae@data, by = "larvae_ID", all = T)
#Remove NAs for when settled and released don't line up
MPA_df <- MPA_df[complete.cases(MPA_df),]
#Need to convert Polygon ID to numeric to sort properly - try to do this earlier in process????
MPA_df$FID_NSB_MP.x <- as.numeric(MPA_df$FID_NSB_MP.x)
MPA_df$FID_NSB_MP.y <- as.numeric(MPA_df$FID_NSB_MP.y)
MPA_df <- MPA_df[with(MPA_df, order(FID_NSB_MP.x, FID_NSB_MP.y)), ]

proc.time() - MPA.dataframe.time

rm(xy, Release_df, Settle_df)


########################################################################
########################################################################
########################################################################
###[6] Creating connectivity tables - (down) column donates to (across) row

#As connectivity matrices
Con_table <- table(MPA_df$FID_NSB_MP.x, MPA_df$FID_NSB_MP.y)
write.csv(Con_table, paste0("./Connectivity_between_MPA_Sarah/output_keep/Con_table/Con_table", size_reduction_threshold[size_reduction_loop], "year", year[year_time], "_pld", pld[pld_time], ".csv"))

#As dataframe
Con_df <- as.data.frame(Con_table)
Con_df$Var1 <- as.character(Con_df$Var1)
Con_df$Var1 <- as.numeric(Con_df$Var1)
Con_df$Var2 <- as.character(Con_df$Var2)
Con_df$Var2 <- as.numeric(Con_df$Var2)
Con_df <- Con_df[with(Con_df, order(Var1, Var2)), ]

#Creating directory
pld_directory <- paste0("./Connectivity_between_MPA_Sarah/output_keep/Con_df/pld", pld[pld_time])
dir.create(pld_directory)
size_reduction_directory <- paste0(pld_directory, "/size_reduction_threshold", size_reduction_threshold[size_reduction_loop])
dir.create(paste0(size_reduction_directory))

##write out Con_df
write.csv(Con_df, paste0(size_reduction_directory, "/Con_df", size_reduction_threshold[size_reduction_loop],"year", year[year_time], "_pld", pld[pld_time], ".csv"), row.names = F)

} #closing pld loop

print(paste0("Finished year", year[year_time]))
} #closing year loop


########################################################################
########################################################################
########################################################################
###[7] Merging connectivity dataframes across years to get average connectivity over decade
#Loading all connectivity dataframes
for (i in 1:length(pld)){
  
  filenames <- list.files(path= paste0("./Connectivity_between_MPA_Sarah/output_keep/Con_df/pld", pld[i], "/size_reduction_threshold", size_reduction_threshold[size_reduction_loop]), pattern= ".csv", full.names=TRUE, recursive=T)

  # load all files into a list
  datalist <- lapply(filenames, read.csv)
  head(datalist)
  
  # set the names of the items in the list, so that you know which file it came from
  datalist <- setNames(datalist,filenames)
  
  # Merging the dataframe (rbind the list)
  dataset <- data.table::rbindlist(datalist, idcol="filename")
  rm(datalist)
  dataset
  
  #Averaging
  Con_df <- group_by(dataset, Var1, Var2)
  Con_df <- dplyr::summarise(Con_df, mean(Freq), sd(Freq))
  
  Con_df <- plyr::rename(Con_df, c("mean(Freq)" = "mean_Freq", "sd(Freq)" = "sd_Freq"))
  
  #Getting rid of rows in dataframe where there is no connection - these exist because they were needed for display in a table
  Nonzero_Con_df <- Con_df[which(Con_df$mean_Freq > 0),]
  
  #Converting to connectivity matrices
  Con_table_mean <- xtabs(Con_df$mean_Freq ~ Con_df$Var1+Con_df$Var2)
  Con_table_mean <- as.data.frame.matrix(Con_table_mean) #this converts to dataframe but might actually be fine to keep as table
  
  Con_table_sd <- xtabs(Con_df$sd_Freq ~ Con_df$Var1+Con_df$Var2)
  Con_table_sd <- as.data.frame.matrix(Con_table_sd) #this converts to dataframe but might actually be fine to keep as table
  
  #creating directory for results output
  mean_results_directory <- paste0("./Connectivity_between_MPA_Sarah/output_keep/Results/pld", pld[i])
  dir.create(mean_results_directory)
  size_reduction_results_directory <- paste0(mean_results_directory, "/size_reduction_threshold", size_reduction_threshold[size_reduction_loop])
  dir.create(paste0(size_reduction_results_directory))
  
  #Writing out connectivity tables and dataframes for mean and standard deviation
  write.csv(Con_df, paste0(size_reduction_results_directory, "/Con_df", size_reduction_threshold[size_reduction_loop], "_mean_pld", pld[i], ".csv"))
  
  write.csv(Nonzero_Con_df, paste0(size_reduction_results_directory, "/NonzeroCon_df", size_reduction_threshold[size_reduction_loop], "_mean_pld", pld[i], ".csv"))
  
  write.csv(Con_table_mean, paste0(size_reduction_results_directory, "/Con", size_reduction_threshold[size_reduction_loop], "_mean_pld", pld[i], ".csv"))
  
  write.csv(Con_table_sd, paste0(size_reduction_results_directory, "/Con", size_reduction_threshold[size_reduction_loop], "_sd_pld", pld[i], ".csv"))
  
  ### [7b] Connectivity dataframes in terms of percentage
  #still needs work to loop and clean
  Con_table_percent <- Con_table_mean
  Released_each_MPA <- Released_larvae@data
  Released_each_MPA <- plyr::count(Released_each_MPA$FID_NSB_MP)
  Released_each_MPA <- dplyr::rename(Released_each_MPA, FID_NSB_MP = x, Number_larvae_release = freq)

  
  for (j in Released_each_MPA$FID_NSB_MP){
    
    Old_row_value <- Con_table_percent[which(row.names(Con_table_percent) == j),]
    New_row_value <- Old_row_value
    
    Con_table_percent[which(row.names(Con_table_percent) == j),] <- (100*Con_table_percent[which(row.names(Con_table_percent) == j),])/(Released_each_MPA$Number_larvae_release[Released_each_MPA$FID_NSB_MP == j])
    
    write.csv(Con_table_percent, paste0(size_reduction_results_directory, "/Con", size_reduction_threshold[size_reduction_loop], "_percent_pld", pld[i], ".csv"))
    
  } #close 7[b] Connectivity dataframes in terms of percentage
  
} #close [7] Merging connectivity dataframes across years to get average connectivity over decade


} #close size_reduction_threshold loop

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


#EXTRA

#head(Con_table_mean)
#Con_table_mean <- Con_table_mean[-1]


