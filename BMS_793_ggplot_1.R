## load ggplot2 package
library(ggplot2)


## Scatterplot of Sepal.Length vs Sepal.Width
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width)) +   ## create plot  
  geom_point()                                       ## add points

## Scatterplot of Sepal.Length vs Sepal.Width for only species setosa
ggplot(subset(iris,Species=="setosa"), aes(x=Sepal.Length, y=Sepal.Width)) +     
  geom_point()  

## Scatterplot of Sepal.Length vs Sepal.Width color coded by species with trend line.
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width,color=Species)) +   
  geom_point() + 
  geom_smooth(method = "lm")          ## add trend line                               

## Side-by-side boxplots of Petal.Length by species.
ggplot(iris, aes(x=Species, y=Petal.Length)) +     
  geom_boxplot()  

## Side-by-side boxplots of Petal.Length by species, but use the log of Petal.Length.
ggplot(iris, aes(x=Species, y=log(Petal.Length))) +     ## uses natural log
  geom_boxplot()  

## Density plot of Sepal.Length.
ggplot(iris, aes(x=Sepal.Length)) + geom_density()

## Density plots of Sepal.Length for each species (all on the same axis system).
ggplot(iris, aes(x=Sepal.Length, color=Species)) + geom_density()

## Density plots of Sepal.Length for each species (all on the same axis system) with fill.
## Note the fill arguement.  alpha makes the less opaque.  Needs to be a value between 0 and 1.
ggplot(iris, aes(x=Sepal.Length, fill=Species)) + geom_density(alpha=0.5)

## Bubble plot of Sepal.Length vs Sepal.Width color coded by species, where the diameter of the bubbles is Petal.Length.
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width,color=Species, size=Petal.Length)) +   
  geom_point(alpha=0.5)  
  
## Panel four plots in one figure.  
## First load (need to install first) the gridExtra package
library(gridExtra)

## Then create the plots.  Note that I'm saving them to variables
p1 <- ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width,color=Species)) +   
  geom_point() + 
  geom_smooth(method = "lm")     

p2 <- ggplot(iris, aes(x=Species, y=Petal.Length)) +     
  geom_boxplot()  

p3 <- ggplot(iris, aes(x=Sepal.Length, fill=Species)) + geom_density(alpha=0.5)

p4 <- ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width,color=Species, size=Petal.Length)) +   
  geom_point(alpha=0.5)  

## use the grid.arrange() function in the gridExtra package
grid.arrange(p1, p2, p3, p4, nrow=2)
## Not pretty, but we can tweak it later. 




