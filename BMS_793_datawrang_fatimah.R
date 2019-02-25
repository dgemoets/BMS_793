library(ggplot2);library(dplyr)

fatim.data <- readRDS("/Volumes/New Volume/flash_2.8.2019/BMS_793B/data.2.13.2019/fatimah.data.rds")

## plot L/R eye
ggplot(fatim.data, aes(x=time, y=amp, color=eye))+geom_point()

## mean within each eye over all time
fatim.data %>% group_by(mouse.id,eye) %>% summarise(mean=mean(amp))

## mean between L/R eye for each time point
fatim.data %>% group_by(mouse.id,time) %>% summarise(mean=mean(amp))

## subset of time >= 0
fatim.data.postme <- subset(fatim.data, time > -1)

## shift amp by amp at t=0
fatim.data.postme$shft.amp <- NA
fatim.data.postme %>% group_by(id) %>% mutate(shft.amp = amp-amp[time==0]) -> tempp

## find max and min (do we need eye?)
fatim.data %>% group_by(id) %>% summarise(max=max(amp))
fatim.data %>% group_by(id) %>% summarise(min=min(amp))

#################################################
## DGs data
dg.data <- read.csv("/Volumes/New Volume/flash_2.8.2019/BMS_793B/data/BMS_793_dg_abund_data.csv")

## separate tax column into two
dg.data$tax2 <- dg.data$tax1 <- NA
temp <- strsplit(as.character(dg.data$tax),"[;]")

for (i in 1:5){
  dg.data$tax1[i] <- temp[[i]][1]
  dg.data$tax2[i] <- temp[[i]][2]
}

temp <- colSums(dg.data[,-c(1,6,7)])  ## 1,6,7 have the taxa  
dg.data.prop <- dg.data  # defining a data frame
dg.data.prop[,-c(1,6,7)] <- sweep(dg.data[,-c(1,6,7)],MARGIN=2,temp,"/")

dg.data.tax1 <- aggregate(.~tax1,dg.data.prop[,-c(1,6)],sum)
dg.data.tax2 <- aggregate(.~tax2,dg.data.prop[,-c(1,7)],sum)

  