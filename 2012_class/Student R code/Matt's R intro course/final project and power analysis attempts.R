setwd("/Users/Hillabeans/Desktop/Academics baby/2009 Fall/R course/for GxE")
library(gregmisc)
data<-read.xls("PlantsAnimals.xls")
mean(data[,2])
sd(data[,2])
rnorm(20,7.85,2)
hist(rnorm(100,7.85,2))
hist(data[,4])

hyb<- data[,6]==1
hyb
means<-aggregate(data, by=list(data[,6]), FUN=mean)
means

Test=rep(NA, 100)
hyb=rep(NA, 20)
non=rep(NA, 20)
for (n in 1:100) {
	hyb=rnorm(20,9,2)
	non=rnorm(20,7,2)
	Test[n]=t.test(hyb,non)
	}
t.test(hyb, non, p.value)
t.test()
warnings()	

library(vegan)
library(ecodist)
#logistic growth equation:
# x(t+1)= x(t) + r*x(t)
x=1
r=x^(1/2)
z=x+x*r
pop=numeric()
for (i in 1:100) {pop[i]=z[i]}
pop

x=1
z=x+x^2
y=0
for (i in 1:100) {y[i]=z[i]}
y

z=c(1:100)
r=z^2
x=z+z*r
y=0
for (i in 1:100) {y[i]=z[i]+z[i]*z[i]^2}
y

plot(y,x)
