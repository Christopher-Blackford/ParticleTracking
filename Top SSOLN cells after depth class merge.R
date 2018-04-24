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
source("K:/Christopher_PhD/Github/ParticleTracking/Particle_Tracking_subcode/2_Setting_up_hexagon_study_extent.R")

Top_percent <- 0.1
Baseline_pld <- 22
Depth_class <- "Intertidal"

########################################################################
########################################################################
########################################################################
########################################################################
###[1] Get overlaps
pld <- c(Baseline_pld)


for (i in 1:length(pld)){
Marxan_output <- read.table(paste0("./Marxan/From_R/output/", Depth_class, "/pld", pld[i], "/Mar_out_ssoln.dat"))
Marxan_output <- dplyr::rename(Marxan_output, Poly_ID = V1, SSOLN = V2)
Marxan_output$pld <- pld[i]
Marxan_output$included <- 1
Marxan_output$RAND <- sample(1:nrow(Marxan_output), nrow(Marxan_output), replace = FALSE)
Marxan_output <- Marxan_output[with(Marxan_output, order(-SSOLN, RAND)),]
Marxan_output <- Marxan_output[Marxan_output$SSOLN >= Marxan_output[ceiling(Top_percent*nrow(Marxan_output)),]$SSOLN,]

Marxan_output_sp <- sp::merge(ConPoly, Marxan_output, by = "Poly_ID", all.x = TRUE)

assign(paste0(Depth_class, "_pld", pld[i]), Marxan_output)

writeOGR(Marxan_output_sp, dsn = "K:/Christopher_PhD/temp", layer = paste0("Marx_SSOLN", Top_percent, "_pld", pld[i]),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)
}


Marxan_output_df <- lapply(ls(pattern=paste0(Depth_class , "_pld")), function(x) get(x))

#Remove baseline pld info so you can compare
Marxan_output_df[[2]] <- NULL
pld <- pld[pld != Baseline_pld]

########################################################################
########################################################################
###Top 10% crossover between different plds
pld_compare <- NULL

for (i in (1:length(Marxan_output_df))){

pld_compare_loop <- merge(Intertidal_pld22, Marxan_output_df[[i]], by = "Poly_ID")
pld_compare_loop <- 100*nrow(pld_compare_loop)/nrow(Intertidal_pld22)

pld_compare <- append(pld_compare, pld_compare_loop)
}

df <- data.frame(pld, pld_compare)
pld_base_add <- data.frame(Baseline_pld, 100)
names(pld_base_add)<-c("pld","pld_compare")
df <- rbind(df, pld_base_add)
df <- df[with(df, order(pld)),]

plot(df$pld, df$pld_compare, xlab = "PLD", ylab = "Percent similar (compared to pld 22)")
lines(df$pld, df$pld_compare)
#title("Percent of cells shared between similar plds in Marxan solution")

