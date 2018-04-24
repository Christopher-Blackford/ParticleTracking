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

rm(sp.na.omit)

###################Loading required packages:
require(plyr)
require(data.table)
require(tidyverse)
require(rgdal)
require(rgeos)
require(maptools)
require(spatialEco)


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

#Loading my MPA shapefile to get proper projection
MPA_mask <- readOGR("K:/Christopher_PhD/CPAWS/Cleaned_standardized/All_PAs", "MPAS_merged")
My_BC_projection <- MPA_mask@proj4string
row.names(MPA_mask@data) <- MPA_mask@data$CB_ID #Change row.names because FID starts at 0 and you want it to start at 1
head(MPA_mask@data)

#Loading Remi's grid where larvae were released
grid <- readOGR("./cuke_present/StudyExtent/Starting_grid", "grid")
NAD_projection <- proj4string(grid)
proj4string(grid)
#Dissolve into one polygon since so you can change grid dimensions
grid <- spTransform(grid, My_BC_projection) #For some reason not "identical" to My_BC_projection, check later
grid <- gUnaryUnion(grid)

#Intersecting - don't know why this works and ConPoly2 <- grid[Ecozone_mask,] doesn't
MPA_mask <- gBuffer(MPA_mask, byid=TRUE, width=0) #Need to do this to avoid ring intersection
row.names(MPA_mask@data) <- MPA_mask@data$CB_ID #Change row.names because FID starts at 0 and you want it to start at 1
MPA_mask_id <- as.character(MPA_mask@data$CB_ID)
ConPoly <- gIntersection(grid, MPA_mask, byid = TRUE, id = MPA_mask_id) 


#Adding dataframe so you can create a shapefile of new study extent
#Clipped dataframe
ConPoly_ID <- row.names(ConPoly)
ConPoly_ID <- as.numeric(ConPoly_ID)
ConPoly_ID <- as.data.frame(ConPoly_ID)
row.names(ConPoly_ID) <- ConPoly_ID$ConPoly_ID

#Original dataframe
ConPoly_data <- as.data.frame(MPA_mask[ConPoly_ID$ConPoly_ID, ])

MPAS <- SpatialPolygonsDataFrame(ConPoly, ConPoly_ID)
MPAS@data <- plyr::rename(MPAS@data, c("ConPoly_ID" = "CB_ID"))

rm(grid, ConPoly_ID, ConPoly, ConPoly_data, MPA_mask_id)


########################################################################
###Removing some of the MPAs that got sliced to thin clipping to Remi's extent

size_reduction_threshold <- 999 #You are only left with 9ish if you go down to 0

#Compare how much smaller clipped MPA layer is to MPA_mask file
MPA_mask@data$Merged_area <- gArea(MPA_mask, byid = TRUE)

MPA_clipped_size <- MPAS
MPA_clipped_size@data$Clip_Area <- gArea(MPA_clipped_size, byid = TRUE)
MPA_clipped_size <- sp::merge(MPA_clipped_size, MPA_mask@data, by = "CB_ID")
MPA_clipped_size@data$size_reduction <- 100*(1 - MPA_clipped_size@data$Clip_Area/MPA_clipped_size@data$Merged_area)
row.names(MPA_clipped_size@data) <- MPA_clipped_size@data$CB_ID #Change row.names 

Size_reduction_df <- MPA_clipped_size@data

#write.csv(Size_reduction_df, "./Connectivity_between_MPA_Stephen/output_keep/size_reduction.csv", row.names = F)

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
MPAS <- MPA_clipped_size[MPA_clipped_size@data$size_reduction <= size_reduction_threshold,]
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

#year_time <- 1
#pld_time <- 1
#i <- 1

###To Loop: Remove hastags below and add closing brackets, delete temporary year_time and pld_time, and concatate year and pld vectors above
###
for (year_time in 1:length(year)){
  
  for (pld_time in 1:length(pld)){

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

#Add larvae IDs to dataset
Con_df <- dataset
Con_df <- subset(Con_df, select = c(long0, lat0, Z0, long, lat, Z, year, rday))
Con_df$larvae_ID <- row.names(Con_df)

#Now you can remove some large files but only if you want to!
rm(dataset)


########################################################################
########################################################################
########################################################################
########################################################################
###[5] Creating dataframe describing release and settlement of each particle
#Clipping to CB_ID to do points in poly
MPAS <- MPAS_loop[,"CB_ID"]

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
MPA_df$CB_ID.x <- as.numeric(MPA_df$CB_ID.x)
MPA_df$CB_ID.y <- as.numeric(MPA_df$CB_ID.y)
MPA_df <- MPA_df[with(MPA_df, order(CB_ID.x, CB_ID.y)), ]

proc.time() - MPA.dataframe.time

rm(xy, Release_df, Settle_df)


########################################################################
########################################################################
########################################################################
###[6] Creating connectivity tables - (down) column donates to (across) row

#As connectivity matrices
Con_table <- table(MPA_df$CB_ID.x, MPA_df$CB_ID.y)
write.csv(Con_table, paste0("./Connectivity_between_MPA_Stephen/output_keep/Con_table/Con_table", size_reduction_threshold,"year", year[year_time], "_pld", pld[pld_time], ".csv"))

#As dataframe
Con_df <- as.data.frame(Con_table)
Con_df$Var1 <- as.character(Con_df$Var1)
Con_df$Var1 <- as.numeric(Con_df$Var1)
Con_df$Var2 <- as.character(Con_df$Var2)
Con_df$Var2 <- as.numeric(Con_df$Var2)
Con_df <- Con_df[with(Con_df, order(Var1, Var2)), ]

#write out Con_df
write.csv(Con_df, paste0("./Connectivity_between_MPA_Stephen/output_keep/Con_df/Con_df", size_reduction_threshold,"year", year[year_time], "_pld", pld[pld_time], ".csv"), row.names = F)

} #closing pld loop

print(paste0("Finished year ", year[year_time]))
} #closing year loop


########################################################################
########################################################################
########################################################################
###[7a] Merging connectivity dataframes across years to get average connectivity over decade
#Loading all connectivity dataframes

#temporary for rockfish project
year <- c(1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007)
pld <- c(30, 60, 120)

size_reduction_threshold <- 999

for (i in 1:length(pld)){
  filenames <- list.files(path="./Connectivity_between_MPA_Stephen/output_keep/Con_df", pattern=glob2rx(paste0("Con_df", size_reduction_threshold, "year*_pld", pld[i], ".csv")), full.names=TRUE,recursive=T)
  
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
  
  Con_df <- plyr::rename(Con_df, c("Var1" = "MPA_Release", "Var2" = "MPA_Settle", "mean(Freq)" = "mean_Freq", "sd(Freq)" = "sd_Freq"))
  
  #Remove rows that don't represent connection between 
  Nonzero_Con_df <- Con_df[which(Con_df$mean_Freq > 0),]
  
  
  ########################################################################
  ###7[b] Optional - Compensating for poor rockfish habitat
  
  Rockfish_habitat_suit <- raster("./Connectivity_between_MPA_Stephen/Rockfish_project/output_rasters/habitat_quality/Rock_habsuit.tif")
  Rockfish_habitat_suit[is.na(Rockfish_habitat_suit[])] <- 0 

  MPA_RockValue <- extract(x=Rockfish_habitat_suit, y=MPAS, fun = mean)
  MPA_RockValue <- data.frame(MPA_RockValue)
  MPA_RockValue$rows <- 1:nrow(MPA_RockValue)
  MPAS@data$rows <- 1:nrow(MPAS@data)
  
  Rock_df <- base::merge(MPAS@data, MPA_RockValue, by = "rows")
  Rock_df <- Rock_df[c("CB_ID", "MPA_RockValue")]

  Habitat_Con_df <- base::merge(Nonzero_Con_df, Rock_df, by.x = "MPA_Release", by.y = "CB_ID", all = T)
  
  Value_to_larvae_constant <- 1/max(Rockfish_habitat_suit@data@values)
  
  Habitat_Con_df$mean_Rockfish_adj <- Habitat_Con_df$mean_Freq*Habitat_Con_df$MPA_RockValue*Value_to_larvae_constant
  
  
  dir.create(paste0("./Connectivity_between_MPA_Stephen/output_keep/Results/pld", pld[i]))
  dir.create(paste0("./Connectivity_between_MPA_Stephen/output_keep/Results/pld", pld[i], "/size_reduction_threshold", size_reduction_threshold))
  
  write.csv(Habitat_Con_df, paste0("./Connectivity_between_MPA_Stephen/output_keep/Results/pld", pld[i], "/size_reduction_threshold", size_reduction_threshold,
                                   "/MPAs_pld", pld[i], ".csv"), row.names = FALSE)
  }


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

########################################

##Converting to connectivity matrices
#Con_table_mean <- xtabs(Con_df$mean_Freq ~ Con_df$MPA_Release+Con_df$MPA_Settle)
#Con_table_mean <- as.data.frame.matrix(Con_table_mean) #this converts to dataframe but might actually be fine to keep as table

#Con_table_sd <- xtabs(Con_df$sd_Freq ~ Con_df$MPA_Release+Con_df$MPA_Settle)
#Con_table_sd <- as.data.frame.matrix(Con_table_sd) #this converts to dataframe but might actually be fine to keep as table

#creating directory for results output
#mean_results_directory <- paste0("./Connectivity_between_MPA_Stephen/output_keep/Results/pld", pld[i])
#dir.create(mean_results_directory)
#size_reduction_results_directory <- paste0(mean_results_directory, "/size_reduction_threshold", size_reduction_threshold)
#dir.create(paste0(size_reduction_results_directory))

#Writing out connectivity tables and dataframes for mean and standard deviation
#write.csv(Con_df, paste0(size_reduction_results_directory, "/Con_df", size_reduction_threshold, "_mean_pld", pld[i], ".csv"))

#write.csv(Con_table_mean, paste0(size_reduction_results_directory, "/Con", size_reduction_threshold, "_mean_pld", pld[i], ".csv"))

#write.csv(Con_table_sd, paste0(size_reduction_results_directory, "/Con", size_reduction_threshold, "_sd_pld", pld[i], ".csv"))

#write.csv(Habitat_handicapped, paste0(size_reduction_results_directory, "/RockfishHab_df", size_reduction_threshold, "_mean_pld", pld[i], ".csv"))



########################################
### [7c] Connectivity dataframes in terms of percentage
#still needs work to loop and clean

#Con_table_percent <- Con_df
#Released_each_MPA <- Released_larvae@data
#Released_each_MPA <- plyr::count(Released_each_MPA$CB_ID)
#Released_each_MPA <- dplyr::rename(Released_each_MPA, CB_ID = x, Number_larvae_release = freq)


#for (j in Released_each_MPA$CB_ID) {
  
#  Old_row_value <- Con_table_percent[which(Con_table_percent$X == j),]
#  New_row_value <- Old_row_value
  
  ##New row = 100 * old row value / total larvae released
#  Con_table_percent[which(Con_table_percent$X == j),] <- (100*Con_table_percent[which(Con_table_percent$X == j),])/(Released_each_MPA$Number_larvae_release[Released_each_MPA$CB_ID == j])
  
#  write.csv(Con_table_percent, paste0(size_reduction_results_directory, "/Con", size_reduction_threshold, "_percent_pld", pld[i], ".csv"))
  
  ########################################
  
