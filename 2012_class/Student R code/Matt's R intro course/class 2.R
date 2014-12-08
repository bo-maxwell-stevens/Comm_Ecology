# changing the working directory: setwd("") and drag your file into the quotes. check it with getwd()
#dir() lists all items in the wd that you set
#read.csv("data.txt")
#name this function ie. data<-read.csv("data.txt")
colnames(data)
#you tell R where column names are located, but default setting for read.csv is first row in spreadsheet by:
#header=false will bring in the frist row as data and not column names
#matt always manipulates data in R, not the original data bc he doesnt want to tamper with it
data[1:3,1:3]
# can also use read.table - more flexible data raeading file fct than read.csv, eg if your data matrix has periods
read.table("data2.txt", sep=" ")
#the sep argument replaces periods with spaces here which structure the data
#header=TRUE says the first row in the data are the column names
#read.xls is a shortcut (frome greg misc package). first need to install his package. the fct changes your excel file into a csv file automatically and reads it into R
install.packages("gregmisc")
library (gregmisc)
read.xls("data3.xls")
read.xls("Invasive data for R.xls")
dir()

# TO select just 1 column in the data:
data[,1]
#there is another way to access this first column using $ and the name of the column
data$env
#another trick: attach(data) which reads in all column in data matrix and names them as new objects by their column name
attach(data)
env
#now just typing env gives all the values in this first column
#detatch fct removes all these new objects created from the attach fct
detach(data)
env
V1 - this isnt found anymore

# important fct in excel is the SORTING fct, we can do this in R
x=c(9,1,4,2,6,7,5)
#to put these in increasing order, run the sort fct
sort(x)
#to put in decreasing order:
sort(x, decreasing=TRUE)
#this isnt sorting x, but creating a new vector x is a new order
#this operation isnt saved. you can save it by saying x=sort(x)

# ORDER - returns a vector of the position of some vector, like a rank, in increasing order
order(x)
# x is 9 1 4 2 6 7 5 
# order(x) returns 2 4 3 7 5 6 1 . this means that the 2nd position in x is the lowest, and going thru, the 1st (eg 1, at the end of the list) position is the highest rank, that is a 9 in x
order(1:7)

#can use order fct to sort martices
y=c(11:20, 1:10)
# this is just some order that you want beforehand
data[order(y),]
# this is ordering your data matrix according to this y order you want. order fct is returning the order of values (posistion of rank orders) of y, which is fed into the row specification of y
# this is essentailly: data[c(11,12,13,14,15...),].
#if you wanted to order data by multiple columns, you nest those in order fct:
# data[order[,2],data[,3],] within data2, whats the order of data3
data[order(V1,V2),]
data[c(4,10,3),]
# this returns the 4th, 10th and 3rd row of data; if you want to order this, can use order fct, put in the brackets
data<-data[order(y),]
data(,V3:V5)
data[,c(1,2,5)]
#this pulls up column 1, 2, and 5
data[c(1,2),c(5)]
#this pulls up the 2 values corresponding to rows 1 and 2 and column 5
data[order(data$V1),]
# this is ordering col V1 in increasing order down the column. this rearranged the row order to make V1 ordered. this is a good sorting mech.
#sort is for a vector or factor and order is for a matrix
inv<-read.xls("Invasive data for R.xls")
inv
inv[order(2),]
inv[order(inv$Invasive.Rank),]
inv$Invasive.Rank
attach(inv)

# SUMMARY - summary of all columns in the data, useful as data quality check...which of my cols are numeric, have missing values (N/A), and can do basic stats, mean, sd, var, correlation, etc
summary(data)
summary(inv)
data$env<-factor(data$env)
# this turned col env to categorical data (character vector). a factor is categorical data
summary(data)
data[,2]
mean(data[,2])
# a shortcut to give mean for each column:
mean(data)
#this gives an error for a col with non numerical data
sd(data)
# this is the standard deviation of the columns
var(data[,2])
# this fct is the variance of eg. col #2
var(data)
# this returns the variance co-variance matrix; not the same as the mean(data). 
cor(data[,1], data[,2])
#returns the correlation statistic for 2 columns or more, eg cor(data) returns correlations for all of the columns. the default is Pearson's R, but can change to Kendall's tau statistis or Spearman's rank row
cor(data, method="spearman")
cor(data, method="kendall")
#to know how many #s or values in a vector use length fct. or if you dont know how many rows in your matrix
length(data[,2])
length(data[3,])
data[5:length(data[,1]),]
# DIM fct shows how many rows and col there are: returns 20 11 which means 20 rows and 11 col
dim(data)
# APPLY fct:
apply(data,MARGIN=1,FUN=sum)
#MARGIN: 1 is rows and 2 is columns, sum adds together all values in rown, but can do this for mean, sd, var, etc. (basic stats)
apply(data, MARGIN=2, FUN=sd)
apply(data, MARGIN=2, FUN=var)
data
data[1,]
ob1<-data[1,-1]
#this removes the first col (env) bc non numerical
ob1
ob1[ob1>0]<-1
# all the values that were greater than zero has now become 1
ob1
# gives 0s and 1s representing species presence/ansence
sum(ob1)# this counts number of species in ob1 (which is the first row without the env column)

# PLOTTING
plot(data$V1)
#plots value of V1 from first to last (in the vector) or top to bottom (in a matrix)
example(plot)

# PAR fct sets graphical parameters for plot
#mfrow for given window, how many plots are going to be in it
par(mfrow=c(1,2))
plot(data$V1)
plot(data$V2)
# this gives both plots in the same window
abline(h=100)
#this draws a line in the V2 plot at 100
plot
