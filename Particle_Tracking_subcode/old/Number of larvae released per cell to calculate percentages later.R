
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Number of larvae released per cell to calculate percentages later
larvae.release.time <- proc.time()

if (year_time == 1 & pld_time == 1){
  #Determining least number of larvae released 
  Released_larvae_df <- Released_larvae@data
  min_release <- count(Released_larvae_df, Poly_ID)
  min_release <- min(min_release$n)
} else {print("Already calculated minimum larval release per cell")
}

Released_larvae_df <- Habitat_classes
Released_larvae_df$RAND <- sample(1:nrow(Released_larvae_df), nrow(Released_larvae_df), replace=F)
Released_larvae_df <- Released_larvae_df[with(Released_larvae_df, order(Poly_ID.x, RAND)), ]

Released_larvae_df$Larv_code <- NA #Need to fill it with something
Released_larvae_df$Larv_code[1] <- 1
row.names(Released_larvae_df) <- 1:nrow(Released_larvae_df)
Released_larvae_df <- Released_larvae_df[1:10000,]

for (release_count in 2:nrow(Released_larvae_df)){
  if (Released_larvae_df$Poly_ID.x[release_count] == Released_larvae_df$Poly_ID.x[release_count-1]){
    Released_larvae_df$Larv_code[release_count] = Released_larvae_df$Larv_code[release_count-1] + 1}
  else {Released_larvae_df$Larv_code[release_count] = 1}
}

Released_larvae_df2 <- Released_larvae_df[(Released_larvae_df$Larv_code < min_release),]

assign(Habitat_classes_names[i], Released_larvae_df2)

proc.time() - larvae.release.time
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}




#Number of larvae released per cell to calculate percentages later
if (year_time == 1 & pld_time == 1){
  #Determining least number of larvae released 
  Released_larvae_df <- Released_larvae@data
  min_release <- count(Released_larvae_df, Poly_ID)
  min_release <- min(min_release$n)
  assign(paste0(Release_classes_names[i], min_release), min_release)
} else {print("Already calculated minimum larval release per cell")
}
