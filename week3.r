commDat=as.data.frame(read.csv("/home/bonobo/Documents/Comm_Eco/commComp_2000.csv",na.strings=c("","NA")))#The na.strings command means that any blank or NA cells will be given a value of NA. We're making it a data.frame so that it differentiates between strings (i.e. words) and numbers

Hymnptr=commDat[commDat$Order=="Hymenoptera",] #This is read as "all columns of commDat where Order is Hymenoptera"
Hymnptr[,1:8]#let's just look at the first eight columns of this subsetted dataset

Hymnptr=subset(commDat,commDat$Order=="Hymenoptera")
Hymnptr[,1:8]

FremDat=commDat[,substr(colnames(commDat),1,3)=="fre"]
#let's breadk this down
colnames(commDat) #the column names of the dataset

substr(colnames(commDat),1,3) #selecting the first three letters of the column names

help(substr) #check out the parameters needed for substr
substr(colnames(commDat),1,3)=="fre" #this returns a TRUE/FALSE record of column names; TRUE if "fre", FALSE if not
head(FremDat)

FremDatTax=cbind(commDat[,1:6],FremDat)
head(FremDatTax)

FremDatHymn=FremDat[commDat$Order=="Hymenoptera",] #this is read as "all columns of FremDat where Order is Hymenoptera"; the latter is indexed froma different object, commDat
FremDatHymn

SR=function(x){
  return(length(x[x!=0])) #this is read as "return the number of observations in a sample that is not zero"
}

commSR=apply(commDat[,-c(1:6)],2,SR) #apply the SR function by columns (1=rows, 2=columns)
#let's break this down
head(commDat[,-c(1:6)]) #we only wanted to calculate species richness for the sample columns; thus, to apply our SR function, we said "apply this function to all columns except the first 6"

commSR #species richness per sample
?apply
treeCat=c(rep(c("Fre","F1","BC","Narr"),each=10))
boxplot(commSR~treeCat,ylab="Species Richness")

summary(aov(commSR~treeCat))#nope!

tiff("/home/bonobo/Documents/Comm_Eco/Species_Richness_boxplot.tif",width=80,height=80,units="mm",res=150)
par(oma=c(0,0,0,0),mar=c(4,4,0,0))
boxplot(commSR~treeCat,xlab="Tree Type",ylab="Species Richness")
dev.off()
?par
library(vegan)

help(specaccum) #let's look at the parameters and outputs for this function
FremAccum=specaccum(t(FremDat),method="random",perm=100)
#Let's break this down
t(FremDat) #specaccum wants the dataset with species as columns and communities as rows; t(FremDat) is the transpose of our Freemont sub-dataset

FremAccum
#This function randomly selects 1,2,...,10 samples 100 times, and calculates the mean and standard deviation of species richness

plot(FremAccum$sites,FremAccum$richness,xlab="Samples",ylab="Richness") #looks like the species pool for Freemont was not fully sampled
#Let's connect the dots
lines(FremAccum$sites,FremAccum$richness)
plot(FremAccum)
AllAccum=specaccum(t(commDat[,-c(1:6)]),method="random",perm=100)
plot(AllAccum)

tiff("G:/LOCKED/R_code/week03/Species_Accum_Curve.tif",width=120,height=80,units="mm",res=150)
par(oma=c(0,0,0,0),mar=c(4,4,0,0))
plot(AllAccum)
dev.off()