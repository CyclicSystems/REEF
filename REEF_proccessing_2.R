library(tidyverse)
library(lunar)
library(lubridate)
library(measurements)
library(sf)


sightings_unfiltered <- read_tsv("sealy_bahamas/sightings.txt")
surveys_unfiltered <- read_tsv("sealy_bahamas/surveys.txt")
geocodes_2 <- read_tsv("sealy_bahamas/Zone4geog.txt")
species_codes_2 <- read_tsv("sealy_bahamas/REEFspecies.txt")

#----------------- lat/long --------------- #

geocodes_2$lat[geocodes_2$lat == "NULL"] <- NA
geocodes_2 <- na.omit(geocodes_2)

geocodes_2$lon[geocodes_2$geogid == 44020020] <- "-72 20.50"
geocodes_2$lon[geocodes_2$geogid == 41050022] <- "-77 09.60"
geocodes_2$lon[geocodes_2$geogid == 41010035] <- "-79 17.45"
geocodes_2$lon[geocodes_2$geogid == 41010069] <- "-79 09.461"
geocodes_2$lon[geocodes_2$geogid == 41010071] <- "-79 11.500"
geocodes_2$lon[geocodes_2$geogid == 41010072] <- "-79 09.675"
geocodes_2$lon[geocodes_2$geogid == 42010066] <- "-77 33.02"
geocodes_2$lon[geocodes_2$geogid == 42090104] <- "-76 05.075"

geocodes_2$lat[geocodes_2$geogid == 41010035] <- "25 30.55"
geocodes_2$lat[geocodes_2$geogid == 41010071] <- "25 22.100"
geocodes_2$lat[geocodes_2$geogid == 41010072] <- "25 08.526"
geocodes_2$lat[geocodes_2$geogid == 42010066] <- "25 01.63"
geocodes_2$lat[geocodes_2$geogid == 42040147] <- "24 51.000"
geocodes_2$lat[geocodes_2$geogid == 42040171] <- "24 53.21"
geocodes_2$lat[geocodes_2$geogid == 42090104] <- "23 45.096"
#fixes a problem with spaces after the decimal in some values

geocodes_2$lat[geocodes_2$geogid == 42040052] <- NA
geocodes_2$lon[geocodes_2$geogid == 42040052] <- NA
#removes datapoint in the middle of Cuba

geocodes_2$lat <- conv_unit(geocodes_2$lat, from = 'deg_dec_min', to = 'dec_deg')
geocodes_2$lon <- conv_unit(geocodes_2$lon, from = 'deg_dec_min', to = 'dec_deg')
#converts to decimal degrees

#--------------------------- data cleanup and joining ----------------- #

names(surveys_unfiltered) <- c("Form", "Exp", "geogid", "Date", "Stemp", "Btemp", "Btime", "Start", "Visibility", "AverageDepth", "MaxDepth", "Current", "Habitat")
bah_data_2 <- left_join(surveys_unfiltered, filter(sightings_unfiltered, Species == 210), By = Form)
bah_data_2$Abundance[is.na(bah_data_2$Abundance)] <- 0
bah_data_2 <- select(bah_data_2, -BTime, -Family, -Species, -Geozone, - MaxDepth)
#Produces a table with one row per survey, and observed abundacnce of Balistes vetula

bah_data_2$Habitat[bah_data_2$Habitat == 0] <- NA
bah_data_2$AverageDepth[bah_data_2$AverageDepth == 0] <- NA
bah_data_2$AverageDepth[bah_data_2$AverageDepth > 10] <- 10
bah_data_2$Current[bah_data_2$Current == 0] <- NA
bah_data_2$Habitat[bah_data_2$Habitat == 0] <- NA
bah_data_2$Habitat[bah_data_2$Habitat == 12] <- NA
bah_data_2$Habitat[bah_data_2$Habitat == 11] <- NA
bah_data_2$Habitat <- as.factor(bah_data_2$Habitat)
bah_data_2$Stemp[bah_data_2$Stemp == 0] <- NA
bah_data_2$Stemp[bah_data_2$Stemp < 60] <- NA
bah_data_2$Btemp[bah_data_2$Btemp == 0] <- NA
bah_data_2$Btemp[bah_data_2$Btemp < 60] <- NA
bah_data_2$Btime[bah_data_2$Btime < 10] <- NA
bah_data_2$Btime[bah_data_2$Btime > 100] <- NA
bah_data_2$Start[bah_data_2$Start == 0] <- NA
bah_data_2$Visibility[bah_data_2$Visibility == 0] <- NA
bah_data_2$AverageDepth <- ordered(as.factor(bah_data_2$AverageDepth))
bah_data_2$Visibility <- ordered(as.factor(bah_data_2$Visibility))
bah_data_2$Current <- ordered(as.factor(bah_data_2$Current))
#some basic data cleanup and casting 

which(is.na(mdy(bah_data_2$Date)))
#5237  7741  7886  9289  9669 10769 10789 10797 11357 11791 11792 are all null

bah_data_2$Date <- mdy(bah_data_2$Date)
bah_data_2$nDate <- decimal_date(bah_data_2$Date)
bah_data_2$nMonth <- bah_data_2$nDate%%1
bah_data_2$moon_phase <- lunar.phase(bah_data_2$Date)
bah_data_2$vetula_breeding <- as.factor(((lunar.phase(bah_data_2$Date, name = 4) == "Full") & ( month(bah_data_2$Date) > 10 | month(bah_data_2$Date) < 4)))

bah_data_2 <- left_join(bah_data_2, geocodes_2, by = "geogid")
bah_data_2$lat <- as.numeric(bah_data_2$lat)
bah_data_2$lon <- as.numeric(bah_data_2$lon)
bah_data_2$k_fold <- sample(1:5, length(bah_data_2$Form), replace = T)
bah_data_2 <- mutate(bah_data_2, Year = as.factor(year(Date)))

#-------------------- funtions -----------------------------#

#converts categories to categorical mean
SFMA_to_num <- function(data){
  Ftable <- tibble(None = sum(data ==0), 
                   Single = sum(data ==1), 
                   Few = sum(data ==2),
                   Many = sum(data ==3),
                   Abundant = sum(data ==4))
  
  if(sum(Ftable[2:4]) != 0){
    AvgF <- sum(Ftable[2:4]*c(2, 4.16, 10)/sum(Ftable[2:4]))
  }
  if(sum(Ftable[3:5]) != 0){
    AvgM <- sum(Ftable[3:5]*c(11, 33.8, 100)/sum(Ftable[3:5]))
  }
  if(sum(Ftable[4:5]) != 0){
    AvgA <- sum(Ftable[4:5]*c(200, 348)/sum(Ftable[4:5]))
  }
  
  score_vector_to_nums <- function(scores){
    nums <- scores
    nums[nums == 2] <- AvgF
    nums[nums == 3] <- AvgM
    nums[nums == 4] <- AvgA
    return(nums)
  }
  
  return(score_vector_to_nums(scores = data))
}

#converts numerical count to SFMA categories
num_to_SFMA <- function(data){
  cats <- data
  cats[cats < .5] <- 0
  cats[cats>=.5 & cats<1.5] <- 1
  cats[cats>1 & cats <= 10] <- 2
  cats[cats>10 & cats <= 100] <- 3
  cats[cats> 100] <- 4
  
  return(cats)
}

CI_from_SFMA<- function(data, include_observational = TRUE){
  Ftable <- tibble(absent = sum(data ==0), 
                   pressent = sum(data > 0), 
                   Abundant = sum(data ==4))
  if(include_observational){
    if((Ftable[3]/sum(Ftable[1:2]) >= .1)){
      CI <- -0.27 + 3.88/((Ftable[2] -1 )^.39)
    } 
    else{
      CI <- .03 + 2.42/((Ftable[2] -1 )^.48)
    }
  }
  else{
    if((Ftable[3]/sum(Ftable[1:2]) >= .1)){
      CI <- -2.16 + 4.22/((Ftable[2] -1 )^.12)
    } 
    else{
      CI <- .06 + 1.28/((Ftable[2] -1 )^.5)
    }
  }
  return(as.numeric(CI))
}

r_square<- function(expected, observed){
  return (1 - sum((expected - observed)^2)/sum((observed-mean(observed))^2))
}

RMSE <- function(expected, observed){
  return (sum((expected-observed)^2)/length(expected))^1/2
}

avg_error <- function(expected, observed){
  return (mean(expected - observed))
}

#------------------------- exporting lat and long for mapping ----------------------#
write_tsv(na.omit(data.frame(geocodes_2$lat, geocodes_2$lon)), path = "used_geozones.txt")
