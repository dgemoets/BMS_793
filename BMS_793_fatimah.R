library(tidyr)
library(splitstackshape)
library(ggplot2)
#####################################################################
## toy example
## unstacked data example
unit1 <- c(2,3,4,5,4)
unit2 <- c(7,9,4,5,3)
unit3 <- c(5,2,3,6,7)
unit4 <- c(9,5,6,7,8)

## create "unstacked" data frame
unstck.data <- data.frame(unit1,unit2,unit3,unit4)
## gather will stack it
stack.data <- gather(unstck.data, key="unit", value="response")

## combining data frames.
## follows order of df above (unit1,unit2,unit3,unit4)
sex<-c("F","F","M","M")
trt <- c("ctrl","A","ctrl","A")
meta.data <- data.frame(sex,trt)

## add columns for sex and trt.
stack.data$sex <- NA
stack.data[stack.data$unit=="unit1",]$sex <- "F" 


#####################################################################
## Fatimah's data

mouse.data <- read.csv("~/Desktop/flash/BMS_793B/fatimah/BMS_793_fatimah_data.csv")
mouse.data1 <- mouse.data[,seq(3,37,by=2)]                  # extract AMP data
mouse.data2 <- gather(mouse.data1, key="id", value="amp")   # stack AMP data

mouse.meta <- read.csv("~/Desktop/flash/BMS_793B/fatimah/meta_data1.csv", header=FALSE)      # load in meta data file
mouse.meta1 <- mouse.meta[c(1,2,11),seq(2,36,by=2)]         # extract sex, left/right and other useful data
mouse.meta2 <- as.data.frame(t(mouse.meta1))                # transpose (using t()) to put into stacked form

mouse.meta2$sex <- c(rep("M",4),"F","F","M","M",rep("F",10))   # tack on sex variable
colnames(mouse.meta2) <- c("ID","mouse.id","eye","sex")        # rename columns with meaningful values

mouse.meta3 <- expandRows(mouse.meta2, 512, count.is.col = FALSE, drop = TRUE)  # Expand rows to match mouse.data2 dataframe
rownames(mouse.meta3) <- c()                                                    # remove row names
time <- rep(mouse.data[,2],18)                                                  # create vector of time values (repeated)
mouse.data3 <- cbind(mouse.data2, mouse.meta3,time)                             # bind all dataframes together

## test plot
ggplot(mouse.data3, aes(x=time,y=amp,color=mouse.id))+geom_point(alpha=0.2)
fatimah.data <- mouse.data3
saveRDS(fatimah.data,"fatimah.data.rds")  # this will save the final data frame


