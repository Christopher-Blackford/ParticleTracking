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

Depth_class <- "Offshore"
#Average_pld <- 22



Contraction_length <- 25
Top_percent <- 0.1

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Getting all PLDs for a Depth class

###Defining present PLDs of species in depth class
Depth_class_Present_directory <- paste0("./Marxan/From_R/output/hexagon/", Depth_class)
PLDs <- list.files(path=Depth_class_Present_directory)
PLDs <- as.numeric(gsub(pattern = "pld", replacement = "", PLDs))
#Do you want to remove mean PLD from analysis?
#PLDs <- PLDs[-which(PLDs==Average_pld)] 

Depth_class_Present_directory <- paste0("./Marxan/From_R/output/hexagon/", Depth_class, "/pld", PLDs)
Depth_class_Present_PLDs <- list.files(path=Depth_class_Present_directory, pattern = "ssoln", full.names=TRUE)
Depth_class_Present_PLDs


###Applying a contracted PLD in the future
PLDs <- sort(PLDs)
Future_PLDs <- PLDs - Contraction_length

#If contraction leads to a negative PLD, just put the PLD = 1
for (i in 1:length(Future_PLDs)){
  if (Future_PLDs[i] < 0){
    Future_PLDs[i] = 1
    }else{}
}

#Directory management
Depth_class_Future_directory <- paste0("./Marxan_future_AllPLD/From_R/output/hexagon/", Depth_class, "/pld", Future_PLDs)

Depth_class_Future_PLDs <- list.files(path=Depth_class_Future_directory, pattern = "ssoln", full.names=TRUE)
Depth_class_Future_PLDs

#Depth class map
Depth_class_map <- readOGR(paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/", Depth_class), Depth_class)
Depth_class_extent <- Depth_class_map@data

########################################################################
########################################################################
########################################################################
########################################################################
###[2] Determining overlaps between present and future plds

df <- matrix(NA, nrow = length(Depth_class_Present_PLDs), ncol = 1)
df <- as.data.frame(df)

for (i in 1:length(Depth_class_Present_PLDs)){
  
  #current
  File_1 <- read.table(Depth_class_Present_PLDs[i])
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
  
  #future  
  File_2 <- read.table(Depth_class_Future_PLDs[i])
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
  df[i,1] <- df_entry
  }
df$V1 <- df$V1*100

colnames(df) <- paste0("Future PLD contracted ", Contraction_length, " days")
rownames(df) <- paste0("PLD ", PLDs)

write.csv(df, paste0("./Marxan_future_AllPLD/Sensitivity/Present_FutureContracted/", Depth_class, "/", Top_percent, "Sensitivity_Contraction_", Contraction_length, ".csv"))


########################################################################
########################################################################
########################################################################
########################################################################
###[3] Synthesizing various contracted PLDs into a meaninful figure

Current_Future_df <- read.csv(paste0("./Marxan_future/Sensitivity/Present_Future/", Depth_class, "/0.1Sensitivity_vector.csv"))
Current_Future_df <- dplyr::rename(Current_Future_df, PLD = X)
Current_Future_df$PLD <- gsub(pattern = "Current to Future ", replacement = "", x = Current_Future_df$PLD)

Contracted_PLDs <- list.files(path=paste0("./Marxan_future_AllPLD/Sensitivity/Present_FutureContracted/", Depth_class),
                              pattern = "Sensitivity_Contraction_", full.names = TRUE)


for (i in 1:length(Contracted_PLDs)){
  
  Contracted_loop <- read.csv(Contracted_PLDs[i])
  Contracted_loop <- dplyr::rename(Contracted_loop, PLD = X)
  
  Current_Future_df <- merge(Current_Future_df, Contracted_loop, by = "PLD")
  
  }

write.csv(Current_Future_df, paste0("./Marxan_future_AllPLD/Sensitivity/Present_FutureContracted/", Depth_class, "/", Top_percent, "Sensitivity_ContractionTotal.csv"), row.names = FALSE)

#####
####
###
##
#END