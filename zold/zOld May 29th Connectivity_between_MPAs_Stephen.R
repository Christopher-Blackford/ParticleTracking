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
require(data.table)
require(tidyverse)
require(readr)
require(rgdal)
require(rgeos)
require(maptools)
require(plyr)
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
### [2] Choosing year of release and pld you are tracking

# List the particle tracking files for that particular year and pld
year <- 1998
pld <- 30

###To Loop: Remove hastags below and add closing brackets, delete temporary year_time and pld_time, and concatate year and pld vectors above
###
#for (year_time in 1:length(year)){
  
  #for (pld_time in 1:length(pld)){

#temporary year_time and pld_time for non loops
year_time <- 1
pld_time <- 1

########################################################################
########################################################################
########################################################################
###[3] Identifying settlement locations and linking to release locations

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
rm(dataset, rl)


########################################################################
########################################################################
########################################################################
########################################################################
###[4] Setting up study extent you will be using to clip your larval release points to your BC study extent

NAD_projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#All MPAs in the region (merged file that had attribute table reduced in R and merged in Arc)
MPAS_merged <- readOGR("C:/Christopher_PhD/CPAWS/Cleaned_standardized/All_PAs", "MPAS_merged")
row.names(MPAS_merged@data) <- MPAS_merged@data$CB_ID #Change row.names because FID starts at 0 and you want it to start at 1
head(MPAS_merged@data)


###Load in MPAs clipped to Remi's extent
#Can either do long way
#source('C:/Christopher_PhD/Github/ParticleTracking/Connectivity_between_MPA/sub_code/Clipping Merged MPA shapefile to Remis grid.R')

#Or since it takes a lot of time, can just load in this file:
MPAS <- readOGR("./Connectivity_between_MPA_Stephen/study_extent", "MPA_clipped")
row.names(MPAS@data) <- MPAS@data$CB_ID #Change row.names because FID starts at 0 and you want it to start at 1
head(MPAS@data)



########################################################################
########################################################################
########################################################################
########################################################################
###[5] Removing some of the MPAs that got sliced to thin clipping to Remi's extent

#Compare how much smaller clipped MPA layer is to MPAS_merged file
MPAS_merged@data$Merged_area <- gArea(MPAS_merged, byid = TRUE)

MPA_clipped_size <- MPAS
MPA_clipped_size@data$Clip_Area <- gArea(MPA_clipped_size, byid = TRUE)

MPA_clipped_size <- sp::merge(MPA_clipped_size, MPAS_merged@data, by = "CB_ID")
MPA_clipped_size@data$size_reduction <- 100*(1 - MPA_clipped_size@data$Clip_Area/MPA_clipped_size@data$Merged_area)
row.names(MPA_clipped_size@data) <- MPA_clipped_size@data$CB_ID #Change row.names 

head(MPA_clipped_size@data)


Size_reduction_df <- MPA_clipped_size@data

write.csv(Size_reduction_df, "./Connectivity_between_MPA/output/size_reduction.csv")

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

#Removing MPAs that were too clipped by Remi's extent (based on percent loss)
size_reduction_threshold <- 999
MPAS <- MPA_clipped_size[MPA_clipped_size@data$size_reduction < size_reduction_threshold,]
#Clipping to CB_ID to do points in poly
MPAS <- MPAS[,"CB_ID"]

########################################################################
########################################################################
########################################################################
########################################################################
###[6] Creating dataframe describing release and settlement of each particle
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
###[7] Creating connectivity tables - (down) column donates to (across) row

#As connectivity matrices
Con_table <- table(MPA_df$CB_ID.x, MPA_df$CB_ID.y)
write.csv(Con_table, paste0("./Connectivity_between_MPA/output/Con_table", size_reduction_threshold,"year", year[year_time], "_pld", pld[pld_time], ".csv"))

#As dataframe
Con_df <- as.data.frame(Con_table)
Con_df$Var1 <- as.character(Con_df$Var1)
Con_df$Var1 <- as.numeric(Con_df$Var1)
Con_df$Var2 <- as.character(Con_df$Var2)
Con_df$Var2 <- as.numeric(Con_df$Var2)
Con_df <- Con_df[with(Con_df, order(Var1, Var2)), ]


#Find largest donors of the whole network
Donor_df <- NULL
for (i in dim(Con_table)[1]:1){
  
  temp <- (sum(Con_table[i,]) - Con_table[i])
  Donor_df <- rbind(temp, Donor_df)
}
Donor_df <- as.data.frame(Donor_df)
Donor_df$CB_ID <- as.numeric(unique(MPA_df$CB_ID.x))
Donor_df <- Donor_df[order(Donor_df$CB_ID),]
Donor_df <- as.data.frame(Donor_df)
colnames(Donor_df)[1] <- "Donor_larv"
row.names(Donor_df) <- (1:nrow(Donor_df))

##How many were released in each MPA in total?
Released_per_MPA <- data.frame(Released_larvae@data)
Released_per_MPA <- count(Released_per_MPA$CB_ID)

Released_per_MPA$x <- as.character(Released_per_MPA$x)
Released_per_MPA$x <- as.numeric(Released_per_MPA$x)

Released_per_MPA <-  Released_per_MPA[order(Released_per_MPA$x),]



rm(i, temp, MPA.dataframe.time)
########################################################################
########################################################################
########################################################################
###[8] Creating shapefiles with larval dispersal data represented

MPASout <- MPAS

MPASout <- merge(MPAS, Donor_df, by = "CB_ID", all=T)
MPASout@data$Area <- gArea(MPASout, byid = TRUE)

#Write out layer
writeOGR(MPASout, dsn = "./Connectivity_between_MPA/output/shapefiles", layer = paste0("MPA_size", size_reduction_threshold, "year", year[year_time], "_pld", pld[pld_time]), 
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)






#} #closing pld loop

#} #closing year loop
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

