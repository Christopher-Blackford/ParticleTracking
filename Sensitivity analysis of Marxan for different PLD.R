##############################################################################################################################
##############################################################################################################################
###Analysing Marxan results and finding top 10%

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
setwd("K:/Christopher_PhD/Github/ParticleTracking") #Need this to access subcode properly

###################Initialize run with these important parameters
my_resolution <- 10000 #defines raster cell size and controls for biased larvae release
#source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")

Depth_class <- "Nearshore"

Top_percent <- 0.5
Baseline_pld <- 56
Min_pld <- 3
Max_pld <- 120

Depth_class_map <- readOGR(paste0("K:/Christopher_PhD/CH1_MPA/Displaying_study_region/CH1/Depth_class_maps/Depth_class_10km/hexagon/", Depth_class), Depth_class)

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Get overlaps

pld <- c(Baseline_pld, Baseline_pld-1, Baseline_pld-2, Baseline_pld+1, Baseline_pld+2, Min_pld, Max_pld)
pld <- sort(pld)

for (i in 1:length(pld)){
Marxan_output <- read.table(paste0("./Marxan/From_R/output/hexagon/", Depth_class, "/pld", pld[i], "/Mar_out_ssoln.dat"))
Marxan_output <- dplyr::rename(Marxan_output, Poly_ID = V1, SSOLN = V2)
Marxan_output <- sp::merge(Depth_class_map, Marxan_output, by = "Poly_ID", all.x = TRUE)
Marxan_output <- Marxan_output@data
Marxan_output$pld <- pld[i]
Marxan_output$included <- 1
Marxan_output$RAND <- sample(1:nrow(Marxan_output), nrow(Marxan_output), replace = FALSE)
Marxan_output <- Marxan_output[with(Marxan_output, order(-SSOLN, RAND)),]
Marxan_output <- Marxan_output[c(1:round(Top_percent*nrow(Marxan_output))),]

assign(paste0(Depth_class, "_pld", pld[i]), Marxan_output)

#Marxan_output_sp <- sp::merge(Depth_class_map, Marxan_output, by = "Poly_ID", all.x = TRUE)
#writeOGR(Marxan_output_sp, dsn = "K:/Christopher_PhD/temp", layer = paste0("Marx_top", Top_percent, "_pld", pld[i]),
 #      driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

}

Marxan_output_df <- lapply(ls(pattern=paste0(Depth_class , "_pld")), function(x) get(x))
str(Marxan_output_df)


########################################################################
########################################################################
###Top 10% crossover between different plds
pld_compare <- NULL

for (i in (1:length(Marxan_output_df))){
  
  pld_compare_loop <- merge(get(paste0(Depth_class, "_pld", Baseline_pld)),
                            get(paste0(Depth_class, "_pld", pld[i])), 
                            by = "Poly_ID")
  
  pld_compare_loop <- 100*nrow(pld_compare_loop)/nrow(get(paste0(Depth_class, "_pld", Baseline_pld)))
  
  pld_compare <- append(pld_compare, pld_compare_loop)
}

df <- data.frame(pld, pld_compare)
df <- df[with(df, order(pld)),]

jpeg(paste0("./Marxan/Sensitivity/", Depth_class, "/SSOLN overlap top ", Top_percent*100, " percent.jpg"))
plot(df$pld, df$pld_compare, xlab = "PLD", ylab = paste0("Percent similar compared to pld", Baseline_pld))
lines(df$pld, df$pld_compare)
title(paste0("Top ", Top_percent*100, "% of cells shared between similar plds in Marxan solution"))
dev.off()














