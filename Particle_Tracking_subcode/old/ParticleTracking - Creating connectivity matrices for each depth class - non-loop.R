#[5] Creating connectivity matrices for each depth class


#NON-LOOP



#Creating separate shapefiles for intertidal, nearshore, and offshore
#########################################################################################################################################
#########################################################################################################################################
#Intertidal

Intertidal_ConPoly <- ConPoly
Intertidal_SR <- NULL


for (i in (unique(Intertidal_Con$Poly_ID.x))){
  
  #Number of particles released that ended up in same cell
  Dispersed_settle <- count(Intertidal_Con$Poly_ID.x == i & Intertidal_Con$Poly_ID.y == i)
  Dispersed_settle <- Dispersed_settle[2,2]
  #Changing NAs to 0s so you can divide
  Dispersed_settle[is.na(Dispersed_settle)] <- 0
  
  #Number of particles released per cell in total
  Total_dispersed <- count(Intertidal_Con$Poly_ID.x == i)
  Total_dispersed <- Total_dispersed[2,2]
  
  #Percent that were released that ended up in the same cell
  SR_Value <- Dispersed_settle/Total_dispersed
  
  Poly_ID <- i
  
  SR_Value <- data.frame(Poly_ID, SR_Value)
  
  Intertidal_SR <- rbind(SR_Value, Intertidal_SR)
}

Intertidal_SR <- merge(Intertidal_ConPoly@data, Intertidal_SR, by = "Poly_ID", all = TRUE)


#To keep all polygon cells even if no particles released from
Intertidal_SR$SR_Value[is.na(Intertidal_SR$SR_Value)] <- 0
Intertidal_ConPoly@data <- merge(Intertidal_ConPoly@data, Intertidal_SR, by = "Poly_ID")

#To keep only polygon cells that have particles released from them
Intertidal_ConPoly@data <- merge(Intertidal_ConPoly@data, Intertidal_SR, by = "Poly_ID")
Intertidal_ConPoly <- sp.na.omit(Intertidal_ConPoly, margin=1)

#Write out polygon
Int_directory <- paste0("./output/shapefiles/Int_SR_pld", pld, "_year", year)
dir.create(Int_directory)
writeOGR(Intertidal_ConPoly, dsn = Int_directory, layer = paste0("Int_hab_pld", pld, "_year", year),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

#########################################################################################################################################
#########################################################################################################################################
#Nearshore
Nearshore_ConPoly <- ConPoly
Nearshore_SR <- NULL

for (i in (unique(Nearshore_Con$Poly_ID.x))){
  
  #Number of particles released that ended up in same cell
  Dispersed_settle <- count(Nearshore_Con$Poly_ID.x == i & Nearshore_Con$Poly_ID.y == i)
  Dispersed_settle <- Dispersed_settle[2,2]
  #Changing NAs to 0s so you can divide
  Dispersed_settle[is.na(Dispersed_settle)] <- 0
  
  #Number of particles released per cell in total
  Total_dispersed <- count(Nearshore_Con$Poly_ID.x == i)
  Total_dispersed <- Total_dispersed[2,2]
  
  #Percent that were released that ended up in the same cell
  SR_Value <- Dispersed_settle/Total_dispersed
  
  Poly_ID <- i
  
  SR_Value <- data.frame(Poly_ID, SR_Value)
  
  Nearshore_SR <- rbind(SR_Value, Nearshore_SR)
}

Nearshore_SR <- merge(Nearshore_ConPoly@data, Nearshore_SR, by = "Poly_ID", all = TRUE)

#To keep all polygon cells even if no particles released from
#Neashore_SR$SR_Value[is.na(Neashore_SR$SR_Value)] <- 0
#Nearshore_ConPoly@data <- merge(Nearshore_ConPoly@data, Neashore_SR, by = "Poly_ID")

#To keep only polygon cells that have particles released from them
Nearshore_ConPoly@data <- merge(Nearshore_ConPoly@data, Nearshore_SR, by = "Poly_ID")
Nearshore_ConPoly <- sp.na.omit(Nearshore_ConPoly, margin=1)

#Write out polygon
Near_directory <- paste0("./output/shapefiles/Near_SR_pld", pld, "_year", year)
dir.create(Near_directory)
writeOGR(Nearshore_ConPoly, dsn = Near_directory, layer = paste0("Near_hab_pld", pld, "_year", year),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)


#########################################################################################################################################
#########################################################################################################################################
#Offshore
Offshore_ConPoly <- ConPoly
Offshore_SR <- NULL

for (i in (unique(Offshore_Con$Poly_ID.x))){
  
  #Number of particles released that ended up in same cell
  Dispersed_settle <- count(Offshore_Con$Poly_ID.x == i & Offshore_Con$Poly_ID.y == i)
  Dispersed_settle <- Dispersed_settle[2,2]
  #Changing NAs to 0s so you can divide
  Dispersed_settle[is.na(Dispersed_settle)] <- 0
  
  #Number of particles released per cell in total
  Total_dispersed <- count(Offshore_Con$Poly_ID.x == i)
  Total_dispersed <- Total_dispersed[2,2]
  
  #Percent that were released that ended up in the same cell
  SR_Value <- Dispersed_settle/Total_dispersed
  
  Poly_ID <- i
  
  SR_Value <- data.frame(Poly_ID, SR_Value)
  
  Offshore_SR <- rbind(SR_Value, Offshore_SR)
}

Offshore_SR <- merge(Offshore_ConPoly@data, Offshore_SR, by = "Poly_ID", all = TRUE)

#To keep all polygon cells even if no particles released from
#Offshore_SR$SR_Value[is.na(Offshore_SR$SR_Value)] <- 0
#Offshore_ConPoly@data <- merge(Offshore_ConPoly@data, Offshore_SR, by = "Poly_ID")

#To keep only polygon cells that have particles released from them
Offshore_ConPoly@data <- merge(Offshore_ConPoly@data, Offshore_SR, by = "Poly_ID")
Offshore_ConPoly <- sp.na.omit(Offshore_ConPoly, margin=1)

#Write out polygon
Off_directory <- paste0("./output/shapefiles/Off_SR_pld", pld, "_year", year)
dir.create(Off_directory)
writeOGR(Offshore_ConPoly, dsn = Off_directory, layer = paste0("Off_hab_pld", pld, "_year", year),
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)


