#Edited 2/3/2023 -- New Harbor Seal dataset incorporated from Eric Ward data
# Removed SRKW

library(tidyverse)
library(dplyr)

seal.data<-read.csv('For Analysis_WoodetalData/harborSeals.csv')
seal.data <- seal.data[,c(1,3)]
seal.data <- seal.data %>%
  rename(count = HS.WA.SanJuanIslands) %>%
  filter(!is.na(count))

#loading in data with PCAs that Chelsea and Tim derived
thedata <- readRDS("compiled_data.RDS")
makezeroone <- function(x) min(1,x)
thedata$foo <- sapply(X = thedata$count, FUN = makezeroone)

colnames(thedata)
unique(thedata$para.spc)
thedata<-thedata[which(thedata$para.spc=="anisakis_sp_all"|thedata$para.spc=="contracaecum_sp_all"),]
unique(thedata$fish.spc)
thedata<-thedata[which(thedata$fish.spc=="rockfish"|thedata$fish.spc=="hake"|thedata$fish.spc=="pollock"|thedata$fish.spc=="smelt"|thedata$fish.spc=="herring"),]
all.dat.wide<-pivot_wider(data = thedata,
                          names_from = para.spc,
                          values_from=count)

#have a dataset with all necessary driver columns but marine mammals

####edit harbor Seal dataset - 
#remove unneccessary columns
#harborseals<-harborseals[,c(1:3,5)]



#now we have one marmam dataset with interpolated values for each species' abundance, merged by year.
#now we need to put together all the fish datasets, and then merge them with this marine mammal dataset

#harborseals<-harborseals %>% rename(year=Year)
#next thing to do is merge the AllFish dataset with the marmam2 dataset to get all the drivers in one MEGA dataframe
#thedata_wide<-merge(all.dat.wide, harborseals,by="year", all.x=T)
#thedata<-merge(thedata,harborseals, by="year", all.x=T) #will need this for later functions

#read in seabird data
seabird<-read.csv("BirdData.csv")
seabird<-seabird %>% rename(year="Year")
#thedata_wide<-merge(thedata_wide, seabird, by="year", all.x=T)
#thedata<-merge(thedata,seabird, bt="year", all.x=T)

thedata$slat <- scale(thedata$lat)
thedata$slat[is.na(thedata$slat)] <- 0
thedata$slong <- scale(thedata$long)
thedata$slong[is.na(thedata$slong)] <- 0

#Going to use the datasets that Wood et al derived, because they have relevant driver data attached.
#reads in the other driver data to merge with the fish datasets
#making one MEGA SPREADSHEET
### Setup data and functions ####
back.convert <- function(x, mean, sd) x * sd + mean 


# Filter out all parasites except Contracaecum, and remove rows with no temperature
spc.data <- dplyr::filter(thedata, para.spc == "contracaecum_sp_all", !is.na(temp_na_rm))

### Select host fish species commonly infected ####
spc.data$fish.spc <- as.character(spc.data$fish.spc)
foo.min <- 0.04
mean.host <- thedata %>%
  group_by(fish.spc) %>%
  summarise(ntotal = n(), npresent = sum(foo), meancount = mean(count))
mean.host$pfoo <- mean.host$npresent / mean.host$ntotal

host.2.keep <- mean.host$fish.spc[mean.host$pfoo>=foo.min]

spc.data <- spc.data %>%
  filter(fish.spc %in% host.2.keep)


spc.data$fish.spc <- as.factor(spc.data$fish.spc)

#### Do year indexing ####
min.year <- min(spc.data$year)
max.year <- max(spc.data$year)
yearlist <- min.year:max.year
year.lookup <- function(x,y) which(y == x)
yearindex <- sapply(X = spc.data$year, FUN = year.lookup, y = yearlist)
nyears <- max.year - min.year + 1


syearindex <- sapply(X = seal.data$Year, FUN = year.lookup, y = yearlist)
ns <- nrow(seal.data)

### Get Temperature Data
temp.data <- read.csv("For Analysis_WoodetalData/RR_temp.csv", header = T)
temp.data <- temp.data%>% 
  dplyr::select(YEAR, temp_na_rm) %>%
  rename(year = YEAR) %>%
  filter(!is.na(temp_na_rm), year <=2018, year >1921)

tyearindex <- sapply(X = temp.data$year, FUN = year.lookup, y = yearlist)

#### Create X and U ####
Xij <- model.matrix(~ -1 + slat + slong + length, data = spc.data)
Uij <- model.matrix(~ -1 + fish.spc, data = spc.data)
ndata <- nrow(spc.data)

x1obs <- seal.data
x2obs <- temp.data

colnames(x1obs) <- c("year", "xt")
colnames(x2obs) <- c("year", "xt")
