#point.in.poly

my.point.in.poly <- function (x, y, poly.id = NULL) 
{
  if (!inherits(y, "SpatialPolygonsDataFrame")) 
    stop("y must be a SpatialPolygonsDataFrame object")
  if ((inherits(x, "SpatialPointsDataFrame") | inherits(x, 
                                                        "SpatialPoints")) == FALSE) 
    stop("x must be a SpatialPointsDataFrame object")
  if (!is.null(poly.id)) {
    if (length(unique(y@data[, poly.id])) != nrow(y)) 
      stop("poly.id not unique")
  }
  if (is.null(poly.id)) {
    y@data$pids <- rownames(y@data)
    poly.id = "pids"
  }
  z <- x[!is.na(sp::over(x, sp::geometry(y))), ]
  z@data <- data.frame(sp::over(x, y))
  row.ids <- rownames(z@data)
  z@data <- data.frame(z@data[, which(names(y) %in% poly.id)])
  names(z) <- poly.id
  rownames(z@data) <- row.ids
  z@data <- data.frame(z@data, merge(z, x@data, by = "row.names", 
                                     all.x = FALSE))
  rm.idx <- which(names(z) %in% c("x", "y", "Row.names", "optional", 
                                  paste0(poly.id, ".1")))
  if (sum(rm.idx) > 0) 
    z@data <- z@data[, -rm.idx]
  z@proj4string <- x@proj4string
  return(z)
}
my.point.in.poly

test_output <- my.point.in.poly(Released_larvae, ConPoly)




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
  
  return(z)
}



#OLD

#Confirming output
Released_larvae_test <- my.point.in.poly(Released_larvae, ConPoly)


writeOGR(Released_larvae_test, dsn = "K:/Christopher_PhD/temp", layer = "test_my_function",
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)

x <- Released_larvae
y <- ConPoly
poly.id <- "Poly_ID"

head(x@data)
head(y@data)
head(z@data)
View(z@data)

#Testing function
z <- x[!is.na(sp::over(x, sp::geometry(y))), ] #clips spatial points to study extent? = 866463 released points but some still have NAs
z_df <- data.frame(sp::over(x, y))
z_df <- na.omit(z_df)

z@data <- z_df
x_df <- x@data
z <- sp::merge(z, x_df, by = "row.names", all.x = FALSE) #and this?
rownames(z@data) <- z@data$Row.names
z@data$Row.names <- NULL


plot(z)
