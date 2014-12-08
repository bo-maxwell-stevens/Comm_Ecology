setwd("/Users/Hillabeans/Desktop/Academics baby/2009 Fall/R course")
library(vegan)
library(ecodist)
library(gregmisc)
data<-read.csv("data.txt")
par(mfrow=c(2,5))
data
plot(data[,1], data[,2], xlabel = env)

# exercise problems:
for (i in 2:11) {plot(data[,1], data[,i], xlab="env", ylab=colnames(data)[i])}
# or for (i in 1:(ncol(data)-1))

# another "pop-quiz" get means of all rows and then get the mean of those means for env of 1 and env of 2
env=data[,1]
V=data[,-1]
data<-as.matrix(data)
r.means<-apply(V, MARGIN=1, FUN=mean)
#apply V with respect to rows (1), and doing mean
# or couldve done r.means=numeric() this is an empty vector, then do a loop: 
for (i in 1:nrow(V)) {r.means[i]=mean(V[i,])}
r.means

#this is how you're supposed to do this (matt's way)
mu=numeric()
for (i in 1:2) {mu[i]=mean(r.means[env==i])}
mu


# MONTE CARLO methods for estimating probs by simulating data

# use random processes to generate data that you can generate probabilities from. Bayesian stats are related to probability, originally for gambling...
# rnorm fct = random sample from a normal distn, eg. rnorm(100, mean=10, sd=5) is 100 samples with these parameters. by default the mean=0 and sd=1 if you dont specify them.
# runif fct for a uniform distn. as random as you can get bc every value has an equal prob of being observed. there are 2 range parameters eg. runif(100, min=0, max=1)
# there are many more types of distn eg poisson
# rt fct is a random sample from a t disnt. pt is the prob from a students t distn
# sample fct pulls a RS from a given set of data values, eg your data vector 
y=sample(x,10)
# round fct rounds-off #s in your data. 

# Coin example
coin = c(0,1)
flips=numeric()
mc=1000 
# this is the number of monte carlo steps. you can call it anything, eg. n
for (i in 1:mc) {flips[i]=sample(coin, 1)}
flips
hist(flips)
flips[flips==0] # want to know how many zeros
length(flips[flips==0])/length(flips)


# 6.2: some probablilty problems in the notes
# 1. craps:
die1=c(1:6)
die2=c(1:6)
n=1000 
rolls=numeric()
for (i in 1:n) {rolls[i]=sample(die1,1)+sample(die2,1)}
length(rolls[rolls==7 | rolls==11])/length(rolls)
length(rolls[rolls==11])/length(rolls)
0.179+0.048
#2. roulette
red=c(1,3,5,7,9,)
notred=c(1:20)
wins=numeric()
mc=1000
for (i in 1:mc) {
	wins[i]=sample(red,1)-sample(notred,1)
	}
length(wins)

# applying MC to biological problems (section 6.3)
x=read.csv("/Users/Hillabeans/Desktop/Academics baby/2009 Fall/R course/MCexample")
x=read.csv("/Users/Hillabeans/Desktop/Academics baby/2009 Fall/R course/MCexample")
x
sd(x[,1])
sd(x[,2])
x.sim=c(x[,1],x[,2])
x.sim
mud.sim=numeric()
mc=1000
for (i in 1:mc) {
	x1.sim=sample(x.sim,25) 
	x2.sim=sample(x.sim,25)
	mud.sim[i]=mean(x1.sim)-mean(x2.sim)
	}
mu.obs=mean(x[,1])-mean(x[,2])
hist(mud.sim, xlim=c(-5,5))
abline(v=mu.obs)
length(mud.sim[mud.sim>=mu.obs])/mc

PA<-read.xls("/Users/Hillabeans/Desktop/Academics baby/2009 Fall/R course/PlantsAnimals.xls")
library(vegan)
library(ecodist)
ind=PA[,4:9]
dep=PA[10:12]
dep
ind
dis.ind = vegdist(ind, method = "euc")
dis.dep = vegdist(dep, method = "bray")
mantel(dis.dep~dis.ind, nperm=1000)
#output =      mantelr        pval1        pval2        pval3    llim.2.5% 
#-0.066466044  0.668000000  0.333000000  0.558000000 -0.148808813 
#  ulim.97.5% 
# 0.002025815 