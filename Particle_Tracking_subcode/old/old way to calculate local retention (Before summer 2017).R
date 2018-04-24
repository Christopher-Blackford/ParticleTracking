#Old way to calculate local retention (Before summer 2017)

##########Surface retention
for (z in 1:length(Habitat_list)){
  Habitat_ConPoly <- ConPoly
  Habitat_SR <- NULL
  
  for (i in (unique(as.data.frame(Habitat_list[z])$Poly_ID.x))){
    
    Habitat_df <- as.data.frame(Habitat_list[z])
    #Number of particles released that ended up in same cell
    Dispersed_settle <- count(Habitat_df$Poly_ID.x == i & Habitat_df$Poly_ID.y == i)
    Dispersed_settle <- Dispersed_settle[2,2]
    #Changing NAs to 0s so you can divide
    Dispersed_settle[is.na(Dispersed_settle)] <- 0
    
    #Number of particles released per cell in total
    Total_dispersed <- count(Habitat_df$Poly_ID.x == i)
    Total_dispersed <- Total_dispersed[2,2]
    
    #Percent that were released that ended up in the same cell
    SR_Value <- Dispersed_settle/Total_dispersed
    
    Poly_ID <- i
    
    SR_Value <- data.frame(Poly_ID, SR_Value)
    
    Habitat_SR <- rbind(SR_Value, Habitat_SR)
  }
  
  Habitat_SR <- merge(Habitat_ConPoly@data, Habitat_SR, by = "Poly_ID", all = TRUE)
  
  #To keep all polygon cells even if no particles released from
  #Habitat_SR$SR_Value[is.na(Habitat_SR$SR_Value)] <- 0
  #Habitat_ConPoly@data <- merge(Habitat_ConPoly@data, Habitat_SR, by = "Poly_ID")
  
  #To keep only polygon cells that have particles released from them
  Habitat_ConPoly@data <- merge(Habitat_ConPoly@data, Habitat_SR, by = "Poly_ID")
  Habitat_ConPoly <- sp.na.omit(Habitat_ConPoly, margin=1)
  
  assign(Habitat_list_name[z], Habitat_ConPoly)
}

rm(Dispersed_settle, Total_dispersed, i, Poly_ID, SR_Value, Habitat_df, Habitat_SR)
rm(Habitat_list, Habitat_list_name, Habitat_ConPoly)
