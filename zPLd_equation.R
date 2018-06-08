temperature <- seq(1, 35, by=1)


eq <- exp(3.17) * (temperature/15)^(-1.4 -0.27*log(temperature/15))
plot(temperature, eq)

#Minimum dif
temperature <- 5
CurrentPLD <- exp(3.17) * (temperature/15)^(-1.4 -0.27*log(temperature/15)) #At temperatureerature 5

temperature <- 7
FuturePLD <- exp(3.17) * (temperature/15)^(-1.4 -0.27*log(temperature/15)) #At temperatureerature 7