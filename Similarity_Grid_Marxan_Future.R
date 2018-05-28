##############################################################################################################################
##############################################################################################################################
###Similarity_Grid_Marxan_Future.R
#Code by: Christopher Blackford (christopher.blackford@mail.utoronto.ca)
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

Depth_class <- "Offshore"

Top_percent <- 0.25

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Getting all PLDs for a Depth class

Depth_class_directory <- paste0("./Marxan_Future/From_R/output/hexagon/", Depth_class)

Depth_class_plds <- list.files(path=Depth_class_directory)

Marxan_list <- NULL
for (i in 1:length(Depth_class_plds)){
  Marxan_outputs <- list.files(path=paste0(Depth_class_directory, "/", Depth_class_plds[i]), pattern = "ssoln", full.names=TRUE)  
  Marxan_list <- append(Marxan_outputs, Marxan_list)
  }

Depth_class_map <- readOGR(paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/", Depth_class), Depth_class)
Depth_class_extent <- Depth_class_map@data

########################################################################
########################################################################
########################################################################
########################################################################
###[2] Determining overlaps between all PLDs

df <- matrix(NA, nrow = length(Marxan_list), ncol = length(Marxan_list))
df <- as.data.frame(df)

for (j in 1:length(Marxan_list)){
  
  for (i in 1:length(Marxan_list)){
    File_1 <- read.table(Marxan_list[j])
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
    
    
    
    File_2 <- read.table(Marxan_list[i])
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

colnames(df) <- Depth_class_plds
rownames(df) <- Depth_class_plds

write.csv(df, paste0("./Marxan_Future/Sensitivity/", Depth_class, "/", Top_percent, "Sensitivity_grid.csv"))

