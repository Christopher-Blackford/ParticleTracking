
#####Testing function
my.point.in.poly <- function (x, y){
  
  #Testing function
  z <- x[!is.na(sp::over(x, sp::geometry(y))), ] #clips spatial points to study extent? = 866463 released points but some still have NAs
  z_df <- data.frame(sp::over(x, y))
  z_df <- na.omit(z_df)
  
  z@data <- z_df
  x_df <- x@data
  z <- sp::merge(z, x_df, by = "row.names", all.x = FALSE) #and this?
  rownames(z@data) <- z@data$Row.names
  z@data$Row.names <- NULL
  z@proj4string <- x@proj4string
  rm(z_df, x_df)
  return(z)
}
