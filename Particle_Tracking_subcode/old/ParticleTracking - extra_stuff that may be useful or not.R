#Extra stuff that may be reinserted or not....

#######################################################################################
########################################################################################
########################################################################################
#Long-form creating connectivity dataframes

#########################################################################################################################################
#########################################################################################################################################
#Intertidal release
Intertidal_larvae <- Con_df[Con_df$Z0 >= -15,]
Intertidal_larvae <- Intertidal_larvae[Intertidal_larvae$Z0 < 15,]
#Intertidal settlement
Intertidal_larvae <- Intertidal_larvae[Intertidal_larvae$Z >= -15,]
Intertidal_larvae <- Intertidal_larvae[Intertidal_larvae$Z < 15,]

#Showing where each larvae begings and ends
Release_df <- subset(Intertidal_larvae, select = c(long0, lat0, Z0, larvae_ID))
Settle_df <- subset(Intertidal_larvae, select = c(long, lat, Z, larvae_ID, year, rday))

#Associate released points with where they were released from
xy <- subset(Release_df, select = c(long0, lat0))

Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
Released_larvae <- spTransform(Released_larvae, ConPoly@proj4string) #use your custom BC projection for this


points.in.poly.time <- proc.time()
require(spatialEco)
Released_larvae <- point.in.poly(Released_larvae, ConPoly) #2 minutes
proc.time() - points.in.poly.time


#Associate settled points with where they settled 
xy <- subset(Settle_df, select = c(long, lat))

Settled_larvae <- SpatialPointsDataFrame(coords = xy, data = Settle_df, proj4string = CRS(NAD_projection))
Settled_larvae <- spTransform(Settled_larvae, ConPoly@proj4string) #use your custom BC projection for this

points.in.poly.time <- proc.time()
Settled_larvae <- point.in.poly(Settled_larvae, ConPoly) #2 minutes
proc.time() - points.in.poly.time

#Join dataframes to make precursor to connectivity matrices
Intertidal_Con <- merge(Released_larvae@data, Settled_larvae@data, by = "larvae_ID", all = T)
#Remove NAs for when settled and released don't line up
Intertidal_Con <- Intertidal_Con[complete.cases(Intertidal_Con),]
Intertidal_Con <- Intertidal_Con[with(Intertidal_Con, order(Poly_ID.x, Poly_ID.y)), ]

#########################################################################################################################################
#########################################################################################################################################
#Nearshore release
Nearshore_larvae <- Con_df[Con_df$Z0 > -60,]
Nearshore_larvae <- Nearshore_larvae[Nearshore_larvae$Z0 < -15,] #change to less than or equal to
#Nearshore settlement
Nearshore_larvae <- Nearshore_larvae[Nearshore_larvae$Z > -60,]
Nearshore_larvae <- Nearshore_larvae[Nearshore_larvae$Z < -15,]

#Showing where each larvae begings and ends
Release_df <- subset(Nearshore_larvae, select = c(long0, lat0, Z0, larvae_ID))
Settle_df <- subset(Nearshore_larvae, select = c(long, lat, Z, larvae_ID, year, rday))

#Associate released points with where they were released from
xy <- subset(Release_df, select = c(long0, lat0))

Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
Released_larvae <- spTransform(Released_larvae, ConPoly@proj4string) #use your custom BC projection for this

points.in.poly.time <- proc.time()
Released_larvae <- point.in.poly(Released_larvae, ConPoly) #2 minutes
proc.time() - points.in.poly.time


#Associate settled points with where they settled 
xy <- subset(Settle_df, select = c(long, lat))

Settled_larvae <- SpatialPointsDataFrame(coords = xy, data = Settle_df, proj4string = CRS(NAD_projection))
Settled_larvae <- spTransform(Settled_larvae, ConPoly@proj4string) #use your custom BC projection for this

points.in.poly.time <- proc.time()
Settled_larvae <- point.in.poly(Settled_larvae, ConPoly) #2 minutes
proc.time() - points.in.poly.time

#Join dataframes to make precursor to connectivity matrices
Nearshore_Con <- merge(Released_larvae@data, Settled_larvae@data, by = "larvae_ID", all = T)
#Remove NAs for when settled and released don't line up
Nearshore_Con <- Nearshore_Con[complete.cases(Nearshore_Con),]
Nearshore_Con <- Nearshore_Con[with(Nearshore_Con, order(Poly_ID.x, Poly_ID.y)), ]


#########################################################################################################################################
#########################################################################################################################################
#Offshore release
Offshore_larvae <- Con_df[Con_df$Z0 > -250,]
Offshore_larvae <- Offshore_larvae[Offshore_larvae$Z0 < -60,] #change to less than or equal to
#Offshore settlement
Offshore_larvae <- Offshore_larvae[Offshore_larvae$Z > -250,]
Offshore_larvae <- Offshore_larvae[Offshore_larvae$Z < -60,]

#Showing where each larvae begings and ends
Release_df <- subset(Offshore_larvae, select = c(long0, lat0, Z0, larvae_ID))
Settle_df <- subset(Offshore_larvae, select = c(long, lat, Z, larvae_ID, year, rday))

#Associate released points with where they were released from
xy <- subset(Release_df, select = c(long0, lat0))

Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
Released_larvae <- spTransform(Released_larvae, ConPoly@proj4string) #use your custom BC projection for this

points.in.poly.time <- proc.time()
Released_larvae <- point.in.poly(Released_larvae, ConPoly) #
proc.time() - points.in.poly.time


#Associate settled points with where they settled 
xy <- subset(Settle_df, select = c(long, lat))

Settled_larvae <- SpatialPointsDataFrame(coords = xy, data = Settle_df, proj4string = CRS(NAD_projection))
Settled_larvae <- spTransform(Settled_larvae, ConPoly@proj4string) #use your custom BC projection for this

points.in.poly.time <- proc.time()
Settled_larvae <- point.in.poly(Settled_larvae, ConPoly) #
proc.time() - points.in.poly.time

#Join dataframes to make precursor to connectivity matrices
Offshore_Con <- merge(Released_larvae@data, Settled_larvae@data, by = "larvae_ID", all = T)
#Remove NAs for when settled and released don't line up
Offshore_Con <- Offshore_Con[complete.cases(Offshore_Con),]
Offshore_Con <- Offshore_Con[with(Offshore_Con, order(Poly_ID.x, Poly_ID.y)), ]


######Remove unnecessary variables
rm(xy, Release_df, Settle_df, Settled_larvae, Released_larvae)
rm(Intertidal_larvae, Nearshore_larvae, Offshore_larvae)





#######################################################################################
########################################################################################
########################################################################################
#Connectivity matricies stuff

#Adding polygon matrix IDs to link after in matrix
#Release polygons (polygonx)
Poly.x <- unique(Intertidal_Con$Poly_ID.x)
Poly.x <- data.frame(Poly.x) #make sure POLY_ID is ordered ascending
colnames(Poly.x)[1] <- "Poly.x" #easy way to rename since you only have 1 column
Poly.x$matrixID_Poly.x <- row.names(Poly.x)
Intertidal_Con <- merge(Intertidal_Con, Poly.x, by.x = "Poly_ID.x", by.y = "Poly.x")

#Settled polygons
Poly.y <- unique(Intertidal_Con$Poly_ID.y)
Poly.y <- data.frame(Poly.y) #make sure POLY_ID is ordered ascending
colnames(Poly.y)[1] <- "Poly.y" #easy way to rename since you only have 1 column
Poly.y$matrixID_Poly.y <- row.names(Poly.y)
Intertidal_Con <- merge(Intertidal_Con, Poly.y, by.x = "Poly_ID.y", by.y = "Poly.y")

Intertidal_Con2 <- Intertidal_Con[(Intertidal_Con$Poly_ID.x %in% 1) & (Intertidal_Con$Poly_ID.y %in% 1), ]
nrow(Intertidal_Con2)

Larvae_released_total <- count(Intertidal_Con$matrixID_Poly.x == i)


#Connectivity matricies stuff
########################################################################################
########################################################################################
########################################################################################


#Inshore connectivity matrices
Intertidal_mat_table <- table(factor(Intertidal_Con$Poly_ID.x,levels=min(Intertidal_Con$Poly_ID.x):max(Intertidal_Con$Poly_ID.x)),
                              factor(Intertidal_Con$Poly_ID.y,levels=min(Intertidal_Con$Poly_ID.y):max(Intertidal_Con$Poly_ID.y))
)

Intertidal_table2 <- table(factor(Intertidal_Con$Poly_ID.x), #row = 2955
                           factor(Intertidal_Con$Poly_ID.y)  #col = 630
)
colnames(Intertidal_table2) = Intertidal_Con$Poly_ID.x



write.csv(Intertidal_matrix, "C:/Christopher_MSc/temp/matrix2.csv")

Intertidal_matrix <- as.matrix(Intertidal_mat_table)
df <- as.data.frame(Intertidal_matrix) #This should work for surfact retention.... com back to it




#THIS DOESN'T WORK 1) POLY_ID ISN'T ORDERED CORRECTLY WITH MATRIX_ID; 2) MATRIXID_POLY.X DOESN'T MATCH WITH MATRIXID_POLY.Y
#SO [i,i] doesn't actually correspond to surface retention
###Surface retention
require(plyr)
SR_list <- NULL
for (i in (min(Intertidal_Con$matrixID_Poly.x):max(Intertidal_Con$matrixID_Poly.x))){
  Larvae_settle_in_release <- Intertidal_matrix[i,i]
  
  Larvae_released_total <- count(Intertidal_Con$matrixID_Poly.x == i)
  Larvae_released_total <- Larvae_released_total[2,2]
  
  #temp <- paste0("Intertidal_poly", i) don't even need to assign unique values
  output <- (Larvae_settle_in_release/Larvae_released_total)
  
  #assign(temp, output) don't even need to assign unique values
  #Need some way to combine variables together
  SR_list <- append(SR_list, output)
}

SR_df <- as.data.frame(cbind(SR_list))
colnames(SR_df)[1] <- paste0("SR_", pld) #easy way to rename since you only have 1 column
SR_df$matrix_IDmerge <- row.names(SR_df)

SR_df <- merge(SR_df, Intertidal_Con, by.x = "matrix_IDmerge", by.y = "")

