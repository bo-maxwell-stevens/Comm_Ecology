# ONE SAMPLE t-test
# is the true mean of the pop different from some hypo value (zero). can be unidirectional or bi-directional. 
path.data<-c(0,2,2,3,12,5,2,0,8,22,0,3,6,5,7)
t.test(path.data)
# turns out p value, t statistic, mean, and 95% CI
help(t.test)
t.test(path.data, alternative = "greater")
t.test(path.data, alternative = c("greater"), mu=2)
#assumptions: parent popn normally distributed. to check, plot sample distn. use HISTOGRAM fct
hist(path.data)
#path.data is counting data, so usu get lots of zeros and cant go below zero, therefore you cant have a normal distn. theres another test that gets around this normal disnt assumtion: 

# WILCOXIN SIGNED RANK test (w/ continuity correction)
wilcox.test(path.data)
# if two values or more are the same, they are tied (get a warning message)
help(wilcox.test)

# 2 SAMPLE t-test
sp1.height <-c(11,13,20,10,15,16,25,9,19,14)
sp2.height<-c(22,22,20,14,28,31,24,32,25,28)
hist(sp1.height)
hist(sp2.height)
# assumption of equal variance - look at sd
sd(sp1.height)
sd(sp2.height)
par(mfrow=c(1,2))
#this lets you look at the 2 histograms side by side. can specify other dimensions. 
t.test(sp1.height, sp2.height)
#if paired, write the argument paired = TRUE

# use BARPLOT to display the means
height.means<-c(mean(sp1.height), mean(sp2.height))
barplot(height.means)
height.means
# to change the height of one of the axis: ylim (or xlim)
barplot(height.means, ylim=c(0,30))
# to label your axes: xlab and ylab
barplot(height.means, ylim=c(0,30), xlab="Plant Spp", ylab="Mean Height", names=c("Sp1", "Sp2"))
# names argument to add Sp1 and Sp2 to the x-axis

# ANOVA 
# this structure is a little different than t-test. there are 3 steps: 1. specify model structure = y~x relates x(indep) to y (dependent) 2. feed this into fct 3. need to summarize output into useful stats
setwd("/Users/Hillabeans/Desktop/School!/2009 Fall/R course")
read.xls("data.txt")
read.csv("data.txt")
data<-read.csv("data.txt")
x = factor (data$env)
# factor makes this data vector into categorical data, w/o this R would assume its numberical and ordered
y = data$V1
model = y~x
aov.out = aov(model)
summary(aov.out)
par(mfrow=c(2,2))
plot(aov.out)
# this plot is way more helpful to look at than the summary output. 4 plots: residual v fitted - 50% of data on one side, and 50 on the other. Normal quartile-quartile plot - can test normality assumption - does our data deviate from the 1 to 1 line. 2 other graphs that Lau doesnt know what they mean. 
install.packages("gregmisc")
library("gregmisc")
y.means=c(mean(y[1:10]), mean(y[11:20]))
barplot2(y.means)
# same as regular barplot, but barplot2 allows error bars to be put in very easily. need upper and lower limits for both
barplot2(y.means, ylim=c(0,300), ylab ="y", xlab ="x", names = c("Level 1", "Level 2"))
y.sd = c(sd(y[1:10]), sd(y[11:20]))
eb.upper = y.means + y.sd
eb.lower = y.means - y.sd
# now plug these into your plot by:
barplot2(y.means, ylim=c(0,300), plot.ci = TRUE, ci.u = eb.upper, ci.l = eb.lower, ylab ="y", xlab ="x", names = c("Level 1", "Level 2"))

# LINEAR REGRESSION
summary(reg.out <-lm(V1~V2, data = data))
# this gives residuals, coefficients, intercept, and a t-table
# our data is just called "data" in this case. if it was called something else, it would be the second "data". lm is linear model

#multiple regression
summary(lm(V1~V2+V3, data = data))
#to test for interaction of 2 factors, use an *
summary (lm(V1~V2+V3+V2*V3, data=data))

# ANCOVA - analysis of covariance
summary(ancova.out<- lm(V1~factor(env)+V2, data=data))
# a categorical and a continuous column. showing significance for V2 bc now accounted for variation in env. not the case with our regression model
# to make a plot to show correlation btw V1 and V2
attach(data)
plot(V2~V1)
# this is a bivariate scatter plot. can now plot regression line on top: 
abline(reg.out)
# how to add text to plot: use text fct. the location of the label will be in the coordinates(x,y) that you specify. can also add points w/ points fct
text(x=305, y=240, labels="p=0.16")
# reg.out is kind of like a list, when unlist it spews out tons of data.
unlist(reg.out)
names(reg.out)
#this returns all that is in reg.out in a better organized way. can pull out individual outputs with histogram:
hist(ancova.out$residuals)
#assumption: residuals normal about zero


## COMMUNITY ANALYSIS - multivariate data!
# usu not going to be able to get every spp, so must assess if youve adequetly sampled the sp pool. and then do stats eg effect on composition of community. can do univariate stuff too, eg sp richness, abundance, etc. (namely the apply fct). 
# always run summary fct to check data quality. also useful is pairs fct: plot matrix of correlations of all columns - good to visualize patterns
pairs(data[,1:5])
install.packages("vegan")
install.packages("ecodist")
install.packages("labdsv")
library(vegan)
library(ecodist)
library(labdsv)

# sampling adequacy: Species Accumulation Curves (specaccum fct)
env = data[,1]
com = data[,-1]
# this separates the sp names (env) with the community data: -1 is you want to get rid of the first column
SA=specaccum(com, method="random")
SA
# the output is using random method to select a row of sp abundances in the matrix and pulling out species at random iteratively. can funnel this into the plot fct which is much easier to interpret.
plot(SA)
# now that weve sampled enough, can do stats. eg is V1 and V2 different in terms of sp comp
# first get distance: in multivar space, how far apart are two values
dis.com=vegdist(com)
# default = bray-curtis ( not euclidean)

# MRPP - multi response permutation procedure. non-parametric test. this compares the random stat generated to our observed stat. can have multiple levels in the groups (more than 2) but harder to interpret. 
mrpp(dis.com, group = env, permutations = 1000)
# look at A value: how good are our groupings? +1 = very good, -1 very bad. 0.25 is ecologically relevant. significance of delta is the P-value

#anosim test
anosim(dis.com, grouping =env, permutations=1000)

#PerMANOVA: non para version of MANOVA. the col name after the tilda is what is affecting the data, eg env. can be categorical or continuous. uses monte carlo stats.
adonis(dis.com~env, permutations = 1000)
# can get R-squared values, eg as N increases, community comps get different. can add more relationships, eg + dis.com~___. can also test for interactions among variables.
# Marty J Anderson - lots of info on her website.

