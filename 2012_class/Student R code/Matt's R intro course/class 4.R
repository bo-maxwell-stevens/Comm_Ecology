# to save a file, eg. to your desktop, use WRITE fcts
write.csv(com.cor, file="comcor")

# MANTEL TEST - correlation btw 2 multivariate sets (matricies). use ecodist, not vegan bc has same structure as other tests (tilda). coomon application is correlating geographic or genetic distances with some community distance measure.

library (gregmisc)
read.csv("data.txt")
d <- read.csv("data.txt")
com<-d[,5:11]
com
library (vegan)
library (ecodist)
geo = com[, 1:3]
com = com[,-1:-3]
dis.geo = vegdist(geo, method = "euc")
dis.com = vegdist(com, method = "bray")
# here we did both euclidean distance or bray-curtis distance. the distance metric must have certain properties that work with your data. eg, geo distances use euc. bray-curtis makes non-parametric/community data behave well. now use the mantel fct to correlate these two
mantel(dis.com~dis.geo, nperm=1000)
# the output is: the mantel r first (your measure of corr btw the 2 matricies), the there are 3 p-values: the first is the lower tail (pval1), then the upper tail (pval2), then the non-directional Ho: r=0 (pval3). the last two numbers are the lower and upper CIs for r.  

# Visualizing your data: Heatmaps and Ordinations (heatmap, vegdist, rda, eigen, nmds, nmds.min, legend)

# HEATMAPS - often used to view molecular data (eg microarray). color intensities show abundances of diff species in each observation that have been clustered together based on their distances. it has lots of defaults that we may not want. 
heatmap(as.matrix(com), scale = "none")
# heatmap wont accept non-numeric data, so must call as.matrix. can also say scale = FALSE. this prevents it from centering our data (standardizes data with respect to mean of each row)
#rows = obs, column = species. it shuffles data to get them to cluster well into the dendrograms. eg rows and cols will not be in same order as original data.
?scale
env<-read.csv("data.txt")[,1]
env
heatmap(as.matrix(com), scale= "none", labRow = env)

# ORDINATIONS #1: PCA
d
com <- d[,5:11]
com
t.com = t(com) 
#this flips rows into columns and vice versa to transpose the matrix
pca.out <- prcomp(t.com)
pca.out
# Rotation output is the principal component outputs for each row: PC1-PC7. this is the max number of axes we can have since we have 7 dimensions in our dataset (7 species in com). % variance explained by each axes by taking the eigen value for that axis and dividing by sum total of all eigen values.
names(pca.out)
# refer to rotation with the $ 
pca.axes = pca.out$rotation
plot(pca.axes)
plot(pca.axes, col=env)
# by default this will plot PC1 and PC2. can also change the shape of the points too rather than color.
# shortcoming of pca: relationships among different dimensions are linear in pca. this usu isnt true for community data. used a lot on morphometrics.

# ORDINATIONS #2: NMDS
# more appropriate for community data - much like pca, but doesnt assume linear relationships btw factors that we're using. 
# note: detatch labdsv package before ordinations. should be using ecodist.
library(ecodist)
dis.com=vegdist(com)
dis.com
#default method for vegdist is bray-curtis
nmds.out <- nmds(dmat = dis.com, mindim = 1, maxdim=4, nits = 10)
plot(nmds.out$stress)
# this output (a scree plot) will tell you the min number of dimensions you can have, below the stress cut-off (0.2). eg, where do we get diminishing returns for stress vs. complexity? we want the least number of dimensions for visualizations. then repeat process for the exact dimensionality for more reps:
nmds.out <- nmds(dmat = dis.com, mindim = 2, maxdim=2, nits = 50)
min.out<-nmds.min(nmds.out, d=2)
# this pulls out the configuration w/ lowest dim and you plot it:
plot(min.out, col=as.numeric(plot), pch = as.numeric(plot))
legend(x="bottomleft", legend = c("JH", "JL", "PJL", "PJH", "PL", "PH"), col = c(1:6), pch = c(1:6))

#INDICATOR SPECIES ANALYSIS (ISA) - how faithful is a particular sp to a grouping variable, eg occuring on a different habitat types. usu. used only in community data
library(labdsv)
ISA <-duleg(com, clustering =env)
summary(ISA)
# only brings up sp with sig. p-values and indicator values greater than .25 (metric of group fidelity) which are sig. indicators of clustering
?duleg

#Relativizing - adjust community data, eg. by sp maximum value to weight them equally. not relativing can hide patterns
# eg if 1 sp is really different in abundance
rel.com = decostand(com, method = "max", MARGIN =2)

# DEBUGGING - error messages
# 1. actually read error messages, can look in help files or look online & google it if dont understand error.
# 2. immediately check for typos 
# 3. look for data structure (data type) errors. eg as.numeric for colors - if env was a factor vector instead of a numeric vector it wouldnt work
# 4. if have multiple errors, fix them in order! this follows the flow of the commands
# use debug (and undebug to close it up) fcts - cracks your fct open and runs it line by line eg:
FtoC <- function(F) {return(F-32)*(5/9)}
debug(FtoC)
FtoC(32)
undebug(FtoC)
FtoC

# OPERATORS - like fcts they do things to values, eg +,:,~, /. there are relational operators: binary operators, like < >, that will return TRUE or FALSE. also works for vectors of numbers of = length and matricies of = dims
1>2
# this returns FALSE
a=1
b=2
a==b
a!=b
# == means does it =?   != means not equal
x=1:10
y=1:10
x==y
cbind(x,y)
# this overlays the binary comparison. to pull out all values of x that are 2, insert a relational operator. give me this value or dont give me this value
x==2
cbind(x,x==2)
x[x==2]
# give me all the values in x that are equal to 2
x=c(1,1,1,2,2,2,3,3,3)
x
length(x)
y=c(1,2,1,3,4,5,6,5,4)
y
y[x==1]
mean(y[x==2])
cbind(x,y)
y[x>=2]
cbind(x,y)
x[x==2 | x==3]
# "|" is the "or" operator and  "&" is the and operator
y=c(1,1,1,1,1,2,2,2,2)
y
x[x==2&y==1]
# to find actual positions in the vector where these arguments are, creat a matrix of whole numbers of the length of your vector, x. then do your operation
(1:length(x))
(1:length(x))[x==1]

# more data manipulations with logicals
data<- read.csv("data.txt")
attach(data)
data[env==1,]
detach(data)

# LOOPS! - do fct over and over. can use FOR fct
x=1
for (i in 1:1000) {x=x+1}
x
# this says take x and add 1 to it, then replaces (updates) x with this value and do it again, 999 more times. the final value of x will be 1001.
# the loop creates i, you just have to call it something, eg i. the 1:1000 is a vector that defines the number of times the loop runs and the values of i at every step of the loop. do that for whatever is inside the {}
x=1:10
for (i in 1:length(x)) {x[i]=x[i]+1}
x
# 1:length(x) runs it the length of vector x
# this loops says take the value of x at position i and add 1 to it, then update the value of x at i, and do it again. now x is 2:11
x=1:10
y=0
for (i in 1:length(x)) {y[i] = x[i]+1}
y
# now y is 2:11

