## partial solutions for the sleep study data
## load ggplot2 package
library(ggplot2)
sleep.study <- read.csv("/Volumes/Untitled/BMS_793B/data/SleepStudy.csv")

# 1
sleep.study$gender <- as.factor(sleep.study$Gender)
levels(sleep.study$gender) <- c("female","male")

# 2
aggregate(GPA~LarkOwl,data=sleep.study,mean)
aggregate(GPA~LarkOwl,data=sleep.study,sd)

# 3
ggplot(sleep.study, aes(x=LarkOwl, y=GPA)) +
  geom_boxplot() + geom_point(position=position_jitter(width=0.1)) #+ geom_jitter(position=c(0.1)  

# 4 
ggplot(sleep.study, aes(x=DepressionScore, y=GPA,color=Stress)) +   
  geom_point() + 
  geom_smooth(method = "lm")  + geom_smooth(aes(color=NULL), method="lm")

# 5
sleep.study$stdGPA <- scale(sleep.study$GPA)
(sleep.study$stdGPA < -3) || (sleep.study$stdGPA > 3)

# 6 
ggplot(sleep.study, aes(x=GPA)) +   
  geom_density() + geom_point(aes(x=3+mean(GPA),y=0), shape=8, color="red") + 
  geom_point(aes(x=-3+mean(GPA),y=0), shape=8, color="red") 

# 7
ggplot(sleep.study, aes(x=GPA, fill=LarkOwl)) + geom_density(alpha=0.5)

# 8 
aggregate(ClassesMissed~DepressionStatus,data=sleep.study,sum)

# 9
with(sleep.study, table(ClassesMissed~DepressionStatus+AnxietyStatus))
xtabs(ClassesMissed~DepressionStatus+AnxietyStatus, data=sleep.study)
## check
sum(subset(sleep.study, (DepressionStatus=="normal") & (AnxietyStatus=="normal"))$ClassesMissed)



