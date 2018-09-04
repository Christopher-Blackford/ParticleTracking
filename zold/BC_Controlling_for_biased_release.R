########################################################################
########################################################################
########################################################################
########################################################################
#[4a-b] BC Controlling for biased release

########################################################################  
###[4b] Controlling for biased release
#Need to set upper bound on how many larvae can be released from grid cell since areas closer to shore will have more larvae released
#Don't need to set lower bound since I'm assuming cells where not a lot of larvae released are areas where not a lot of that depth class exists
if (Bias_release_files_preloaded == FALSE & year_time == 1 & pld_time == 1){
  
  Biased_release <- Con_df
  set.seed(1)
  Biased_release$RAND <- sample(1:nrow(Biased_release), nrow(Biased_release), replace=F)
  Biased_release <- Biased_release[with(Biased_release, order(Poly_ID.x, RAND)), ]
  
  Biased_release$Larv_code <- NA #Need to fill it with something
  Biased_release$Larv_code[1] <- 1
  row.names(Biased_release) <- 1:nrow(Biased_release)
  
  ###Counting up number of larvae released in each cell. Takes over a day to do this for all habitat classes
  for (release_count in 2:nrow(Biased_release)){
    if (Biased_release$Poly_ID.x[release_count] == Biased_release$Poly_ID.x[release_count-1]){
      Biased_release$Larv_code[release_count] = Biased_release$Larv_code[release_count-1] + 1}
    else {Biased_release$Larv_code[release_count] = 1}
  }
  
  Biased_release <- Biased_release[c("larvae_ID", "long0", "lat0", "Z0", "Poly_ID.x", "RAND", "Larv_code")]
  write.csv(Biased_release, "K:/Christopher_PhD/Github/ParticleTracking/BC_ConnectivityProject/output/Release_bias/ReleaseBias.csv")
  
  ###when Bias_release_files_preloaded == TRUE
} else{print(paste0("Loading in previous ", Habitat_classes_names[i], " larval release file")) 
  Biased_release <- read.csv("K:/Christopher_PhD/Github/ParticleTracking/BC_ConnectivityProject/output/Release_bias/ReleaseBias.csv")
}

###End of biased release file creation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#In cases where the above loop isn't true, make sure you have release dataframe to merge with habitat dataframe later
###Taking biased release file and un-biasing your depth class files
Released_larvae_df <- merge(Biased_release, Con_df, by = c("larvae_ID", "long0", "lat0", "Z0", "Poly_ID.x")) #this is a dumb way to merge
Released_larvae_df <- Released_larvae_df[(Released_larvae_df$Larv_code <= 6100),] #change this 30 to min_release or something more automated at some point
#Remove NAs for when settled and released don't line up
Released_larvae_df <- Released_larvae_df[complete.cases(Released_larvae_df[,"Poly_ID.y"]),]
Released_larvae_df <- Released_larvae_df[with(Released_larvae_df, order(Poly_ID.x, Poly_ID.y)), ]
#write out final release bias file so it's obvious you only include 100 per cell
write.csv(Released_larvae_df, paste0("K:/Christopher_PhD/Github/ParticleTracking/output/Release_bias/merged_bias_file/Release_bias_controlled.csv"))


###now you have controlled for biased release~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

########################################################################
########################################################################
