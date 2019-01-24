library(ggplot2)
data(iris)

## 1

### c
p.box <- ggplot(iris, aes(x=Species, y=Petal.Length)) + 
  geom_boxplot()
p.box

p.scatter <- ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) + 
  geom_point()
p.scatter

## 4
p.density <- ggplot(iris, aes(x=Sepal.Length, color=Species)) + 
  geom_density()
p.density
