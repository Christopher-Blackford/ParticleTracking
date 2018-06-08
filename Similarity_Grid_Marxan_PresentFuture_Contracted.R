##############################################################################################################################
##############################################################################################################################
###Marxan similarity grid

#
##
###
####
#####
#Clear workspace
rm(list=ls())

require(rgdal)
require(rgeos)
require(tidyverse)

###################Initialize run with these important parameters

Depth_class <- "Intertidal"

Top_percent <- 0.1

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Getting all PLDs for a Depth class

Depth_class_Present_directory <- paste0("./Marxan/From_R/output/hexagon/", Depth_class)
Depth_class_Present_PLDs <- list.files(path=Depth_class_Present_directory)

Depth_class_Future_directory <- paste0("./Marxan_future_AllPLD/From_R/output/hexagon/", Depth_class)
Depth_class_Future_PLDs <- list.files(path=Depth_class_Future_directory)

Marxan_Present <- NULL
for (i in 1:length(Depth_class_Present_PLDs)){
  Marxan_outputs <- list.files(path=paste0(Depth_class_Present_directory, "/", Depth_class_Present_PLDs[i]), pattern = "ssoln", full.names=TRUE)  
  Marxan_Present <- append(Marxan_outputs, Marxan_Present)
  }

Marxan_Future <- NULL
for (i in 1:length(Depth_class_Future_PLDs)){
  Marxan_outputs <- list.files(path=paste0(Depth_class_Future_directory, "/", Depth_class_Future_PLDs[i]), pattern = "ssoln", full.names=TRUE)  
  Marxan_Future <- append(Marxan_outputs, Marxan_Future)
}

Depth_class_map <- readOGR(paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/", Depth_class), Depth_class)
Depth_class_extent <- Depth_class_map@data

########################################################################
########################################################################
########################################################################
########################################################################
###[2] Determining overlaps between all PLDs

df <- matrix(NA, nrow = length(Marxan_Future), ncol = length(Marxan_Present))
df <- as.data.frame(df)

for (j in 1:length(Marxan_Future)){
  
  for (i in 1:length(Marxan_Present)){
    File_1 <- read.table(Marxan_Future[j])
    File_1 <- dplyr::rename(File_1, Poly_ID = V1, SSOLN = V2)
    File_1 <- sp::merge(Depth_class_map, File_1, by = "Poly_ID", all.x = TRUE)
    File_1 <- File_1@data
    File_1$included <- 1
    set.seed(1)
    File_1$RAND <- sample(1:nrow(File_1), nrow(File_1), replace = FALSE)
    File_1 <- File_1[with(File_1, order(-SSOLN, RAND)),]
    File_1 <- File_1[c(1:round(Top_percent*nrow(File_1))),]
    File_1 <- File_1[with(File_1, order(Poly_ID)),]
    File_1 <- File_1[,1]
    
    
    
    File_2 <- read.table(Marxan_Present[i])
    File_2 <- dplyr::rename(File_2, Poly_ID = V1, SSOLN = V2)
    File_2 <- sp::merge(Depth_class_map, File_2, by = "Poly_ID", all.x = TRUE)
    File_2 <- File_2@data
    File_2$included <- 1
    set.seed(1)
    File_2$RAND <- sample(1:nrow(File_2), nrow(File_2), replace = FALSE)
    File_2 <- File_2[with(File_2, order(-SSOLN, RAND)),]
    File_2 <- File_2[c(1:round(Top_percent*nrow(File_2))),]
    File_2 <- File_2[with(File_2, order(Poly_ID)),]
    File_2 <- File_2[,1]
    
    df_entry <- length(intersect(File_1, File_2))/length(File_1)
    df[j, i] <- df_entry
  }
  
}

#df <- as.matrix(df)

colnames(df) <- paste0("Current ", Depth_class_Present_PLDs)
rownames(df) <- paste0("Future ", Depth_class_Future_PLDs)

write.csv(df, paste0("./Marxan_future_AllPLD/Sensitivity/Present_FutureContracted/", Depth_class, "/", Top_percent, "Sensitivity_grid.csv"))

