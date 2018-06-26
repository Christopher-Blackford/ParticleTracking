#Plotting the curve

temperature <- seq(1, 35, by=0.1)


eq <- exp(3.17) * (temperature/15)^(-1.4 -0.27*log(temperature/15))
plot(temperature, eq)

#Intital SST range from 4-12
#Mean temp at bottom range from 3-12 but generally much colder

#SST predicted to increased by 0.5 - 2.0 degrees


#Largest contraction would be at coldest temperature (4 degrees) and warmest shift (+2 degrees)
temperature_current_large <- 4
CurrentPLD <- exp(3.17) * (temperature_current_large/15)^(-1.4 -0.27*log(temperature_current_large/15))

temperature_future_large <- 6
FuturePLD <- exp(3.17) * (temperature_future_large/15)^(-1.4 -0.27*log(temperature_future_large/15)) #At temperatureerature 7

CurrentPLD - FuturePLD #26 day contraction

#Smallest contraction would be at warmest temperature (12 degrees) and coldest shift (+0.5 degrees)

temperature_current_small <- 12
CurrentPLD <- exp(3.17) * (temperature_current_small/15)^(-1.4 -0.27*log(temperature_current_small/15))

temperature_future_small <- 12.7
FuturePLD <- exp(3.17) * (temperature_future_small/15)^(-1.4 -0.27*log(temperature_future_small/15)) #At temperatureerature 7

CurrentPLD - FuturePLD #2 day contraction