setwd("Y:/LOCKED/arthropod_data")#Yours may not be the "U:" drive, make sure you enter the correct letter
commDat=as.data.frame(read.csv("commComp_2000.csv",na.strings=c("","NA")))#The na.strings command means that any blank or NA cells will be given a value of NA. We're making it
head(commDat)
colnames(commDat)=c("ArthroID",colnames(commDat[,-1]))
head(commDat)#you can see the first column now has a name
summary(commDat$Order)
summary(commDat[,2])
summary(commDat$fre1)
summary(commDat$Genus)
length(commDat$Genus[!is.na(commDat$Genus)]=="TRUE")#This tells us the total number of commDat$Genus values that are not NA
length(commDat$Genus[!is.na(commDat$Genus)]=="TRUE")/length(commDat$Genus)#Or the proportion of arthropods keyed to Genus
summary(commDat$Species)
length(commDat$Species[!is.na(commDat$Species)]=="TRUE")
length(commDat$Species[!is.na(commDat$Species)]=="TRUE")/length(commDat$Species)
hist(commDat$fre1,breaks=c(-1:max(commDat$fre1)))
commDat$ArthroID[commDat$fre1==max(commDat$fre1)]
## [1] Blotch mine
totals=rowSums(commDat[,-c(1:6)])
hist(totals,breaks=c(-1:max(totals)))
totals1=apply(commDat[,-c(1:6)],1,sum)#We are applying the function sum() to all rows (1=rows,2=columns) of our dataframe commDat excluding the first six rows; that is how thi
totals1#you can see that this is the same as "totals" above
rowVars=apply(commDat[,-c(1:6)],1,var)
plot(totals,rowVars)
setwd("y:/LOCKED/arthropod_data")
commDat=as.data.frame(read.csv("commComp_2000.csv",na.strings=c("","NA")))
Hymnptr=commDat[commDat$Order=="Hymenoptera",]
Hymnptr[,1:8]
Hymnptr=subset(commDat,commDat$Order=="Hymenoptera")
Hymnptr[,1:8]
FremDat=commDat[,substr(colnames(commDat),1,3)=="fre"]
colnames(commDat)
substr(colnames(commDat),1,3)
help(substr)
substr(colnames(commDat),1,3)=="fre"
head(FremDat)
FremDatTax=cbind(commDat[,1:6],FremDat)
head(FremDatTax)
FremDatHymn=FremDat[commDat$Order=="Hymenoptera",]
FremDatHymn
SR=function(x){
  return(length(x[x!=0]))
}
commSR=apply(commDat[,-c(1:6)],2,SR)
head(commDat[,-c(1:6)])
commSR=apply(commDat[,-c(1:6)],2,SR) #apply the SR function by columns (1=rows, 2=columns)
#let's break this down
head(commDat[,-c(1:6)])
commSR #species richness per sample
treeCat=c(rep(c("Fre","F1","BC","Narr"),each=10))
boxplot(commSR~treeCat,ylab="Species Richness")
summary(aov(commSR~treeCat))#nope!
tiff("y:/OPEN/Prince/Species_Richness_boxplot.tif",width=80,height=80,units="mm",res=150)
par(oma=c(0,0,0,0),mar=c(4,4,0,0))
boxplot(commSR~treeCat,xlab="Tree Type",ylab="Species Richness")
dev.off()
library(vegan)
help(specaccum) #let's look at the parameters and outputs for this function
FremAccum=specaccum(t(FremDat),method="random",perm=100)
t(FremDat) #specaccum wants the dataset with species as columns and communities as rows; t(FremDat) is the transpose of our Freemont sub-dataset
FremAccum
plot(FremAccum$sites,FremAccum$richness,xlab="Samples",ylab="Richness") #looks like the species pool for Freemont was not fully sampled
lines(FremAccum$sites,FremAccum$richness)
plot(FremAccum)
AllAccum=specaccum(t(commDat[,-c(1:6)]),method="random",perm=100)
plot(AllAccum)
tiff("G:/Open/Prince/Species_Accum_Curve.tif",width=120,height=80,units="mm",res=150)
par(oma=c(0,0,0,0),mar=c(4,4,0,0))
plot(AllAccum)
dev.off()
tiff("G:/LOCKED/R_code/week03/Species_Accum_Curve.tif",width=120,height=80,units="mm",res=150)
par(oma=c(0,0,0,0),mar=c(4,4,0,0))
plot(AllAccum)
dev.off()
