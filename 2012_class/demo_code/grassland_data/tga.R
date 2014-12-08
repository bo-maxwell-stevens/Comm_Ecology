# Trait-gradient analysis functions

# Raw data file has 4 or more columns, sep by tabs
#	1-3 labelled 'plot', 'species', 'abund' (in any order)
#	4... plot level trait data, for as many traits as desired
# If trait data at plot level are missing, species means may be
# used from separate file ('sppmeans.txt')

# Option to calculate SMA regression of species slopes uses
# smatr routine - remove comment here to load smatr
#library(smatr)

#load.data reads in data file, adds NA column if traits are 
#missing, names columns and returns dataframe for next step
#trt.col indicates column number to use for trait data (= 4
#if there is only one trait, or higher if file has >1)
load.data = function(datafile,trt.col) {
	d = read.delim(datafile)
	if (is.na(trt.col)) {
			d1 = data.frame(d$species,	
				d$plot,d$abund,rep(NA,nrow(d)))
			names(d1)=c('species','plot','abund','tNA')
		} else {
			d1 = data.frame(d$species,d$plot,d$abund,d[,trt.col])
			names(d1)=c('species','plot','abund'
				,names(d)[trt.col])
		}
	d1 = d1[order(d1$species),]
	return(d1)
	}
	
# function to calculate species means (weighted or unweighted)
# from available datda and replace NA trait data with species 
# means; should only be used if at least one trait measurement is 
# available for every species
expand.spp.means.internal <- function(x,wtd) {
	wtd.sum = tapply(x$abund*x[,4],
		as.factor(x$species),sum,na.rm=TRUE)
	wts = tapply(x$abund[!is.na(x[,4])],
		as.factor(x$species[!is.na(x[,4])]),sum)
	sppmeans = wtd.sum/wts

	#sppmeans = tapply(x[,4],x[,1],mean,na.rm=TRUE)
	spp.missing = sum(is.na(sppmeans))
	if (spp.missing>0) {
		print(paste(spp.missing,' spp has/have no trait data; function aborted'))
		break
		}
	spp.match=match(x[,1],names(sppmeans))
	expand.means=sppmeans[spp.match]
	x[,5]=x[,4]
	names(x)[5] = names(x)[4]
	names(x)[4] = paste(names(x)[4],'.exp',sep='')
	x[is.na(x[,4]),4]=expand.means[is.na(x[,4])]
	return(x)
	}
	
# function to replace NA trait data with species means from 
# separate file; this should be used if there is no trait data in 
# D for one or more species, or if desired species mean trait 
# values are different from means of trait measurements within 
# D. 
expand.spp.means <- function(x,meansfile,colnum) {
	means=read.delim(meansfile)
	spp.match = match(x[,1],means[,1])
	expand.means = means[spp.match,colnum]
	x[,5]=x[,4]
	if (length(x[is.na(x[,4]),4])==nrow(x)) {
		names(x)[5] = names(means)[colnum]
		} else {
		names(x)[5] = names(x)[4]
		}
	names(x)[4] = paste(names(x)[5],'.exp',sep='')
	x[is.na(x[,5]),4]=expand.means[is.na(x[,4])]
	missing = sum(is.na(x[,4]))
	if (missing>0) print(paste('Warning: trait values missing for ',missing,'record(s)')) 
	return(x)	
	}
	
#Calculate site and species means and add them as new columns
#option weighted = TRUE to calculate abundance weighted plot
#and species means; FALSE to weight each species occurrence equally.
#Returns data frame with 7 columns: species, plot, abund, 
#tval (trait data with species means as needed), tp (plot mean
#trait value), ts (species mean trait values), and tvalNA (original
#trait data including NA if not measured locally in each plot).
calc.means <- function(d,weighted=TRUE) {
	
	if (dim(d)[2]==5) expTr=TRUE else expTr=FALSE
	if (!weighted) d$abund = 1

	trt.avail=which(!is.na(d[,4]))
	wtd.sum <- tapply(d[trt.avail,4]*d$abund[trt.avail],as.factor(d$plot[trt.avail]),sum,na.rm=TRUE)
	wts <- tapply(d$abund[trt.avail]
		,as.factor(d$plot[trt.avail]),sum,na.rm=TRUE)
	site.T <- wtd.sum/wts
	site2x <- match(d$plot,names(site.T))

	wtd.sum <- tapply(d[,4]*d$abund,as.factor(d$species),sum)
	wts <- tapply(d$abund,as.factor(d$species),sum)
	spec.T <- wtd.sum/wts
	spec2x <- match(d$species,names(spec.T))
	
	
	T <- cbind(d[,1:4],site.T[site2x],spec.T[spec2x])
	if (expTr) T <- cbind(T,d[,5]) else T <- cbind(T,d[,4])
	names(T)=c('species','plot','abund',
		'tval','tp','ts','tvalNA')
	return(T)
	}
	

#Calculate species attributes
# T has 7 columns listed above
# this function assumes no missing vals in tval
# minsp and maxsp can be set to a subset of taxa if
#	results are not needed for entire dataset; this is used
#	by species jackknife analysis.
# calc.slope = FALSE to turn off intraspecific slopes; this
#	saves times if they are not needed or not available when
#	only species mean trait data is used, and it useful for 
#	bootstrap and null model runs.
# use.spmeans = TRUE to calculate slopes including species means where
#	in situ trait data not available; this is only useful to
#	obtain slopes that will plot neatly with plotTS function;
#	useNA = FALSE to exclude species means and only calculate
#	intraspecific slopes using trait data measured in situ.
# 
# This function returns S dataframe with 22 columns and a row for
#	each species:
#	species: species name
#	NPlots: number of plots species occurs in
#	mean.abund: mean species abundance across plots
#	ts: mean trait value
#	betaT: beta trait value
#	alphaT: alpha trait value
#	Rs: niche breadth
#	slopeN: number of plots in which trait
#		was measured in situ; used for slopes analysis
#	b.inter: intercept of abundance weighted regression
#	bs: slope of weighted regression
#	bs.se: standard error of weighted slope
#	bs.p: significance of weighted slope (relative to 0)
#	u.inter: intercept of unweighted regression
#	us: slope of unweighted regression
#	us.se: standard error of unweighted slope
#	us.p: significance of unweighted slope (relative to 0)
#	sma.slope: slope of standardized major axis of species vs. plot traits
#		sma slope is on unweighted trait values
#	sma.lowci: 2.5% lower CI of sma slope
#	sma.hici: 97.5% upper CI of sma slope
#	min.tp: lowest mean trait value of occupied plots
#	max.tp: highest mean trait value of occupied plots
#	min.tval: lowest measured trait value
#	max.tval: higheset measured trait value
#	low.bound: does species occur in plot at lower edge
#		of gradient (1=yes)
#	hi.bound: does species occur in plot at upper edge

calc.spec.attr <- function(T,minsp=1,
		maxsp=length(unique(T$species)),calc.slopes=TRUE,
		use.spmeans=TRUE) {
	#range of site means trait vals
	plotT.r <- range(T$tp,na.rm=TRUE)
	
	#number of occurrences of species
	sp.N <- tapply(T$ts,as.factor(T$species),length)
	sp.N[which(is.na(sp.N))]=0
	Nsp = length(sp.N)
	
	#mean abundance of species
	sp.abund <- tapply(T$abund,
		as.factor(T$species),mean,na.rm=TRUE)
	
	#site trait mean for sites occupied by species
	#this is not recalculated for bootstrap and null
	#assumes that plot characteristics are fixed
	wtd.sum <- tapply(T$tp*T$abund,
		as.factor(T$species),sum)
	wts <- tapply(T$abund,as.factor(T$species),sum)
	sp.betaT <- wtd.sum/wts
	
	#species trait mean
	#this is already calculated in T but is recalculated
	#here for bootstrap and null models where T may have been
	#resampled
	wtd.sum <- tapply(T$tval*T$abund,
		as.factor(T$species),sum)
	wts <- tapply(T$abund,as.factor(T$species),sum)
	sp.T <- wtd.sum/wts

	#species trait mean
	#sp.T <- tapply(T$ts,
	#	as.factor(T$species),median)
	
	#delta trait value (species mean - occupied site mean)
	sp.alphaT <- sp.T-sp.betaT
	
	sp.Ntval=rep(NA,Nsp)
	sp.NtvalNA=rep(NA,Nsp)
	sp.R=rep(NA,Nsp) #niche breadth
	ncases=rep(NA,Nsp) #number of cases with local trait data
	b=array(rep(NA,Nsp*2)) #weighted intra-slope coefficients
	dim(b)=c(Nsp,2)
	b.se = rep(NA,Nsp) #standard error of weighted slope
	b.p = rep(NA,Nsp) #significance of weighted intra-slope
	u=array(rep(NA,Nsp*2)) #unweighted intra-slope coefficients
	dim(u)=c(Nsp,2)
	u.se = rep(NA,Nsp) #standard error of unweighted slope
	u.p = rep(NA,Nsp) #significance of unweighted intra slope
	sm = matrix(NA,Nsp,3) #smatr slope, lo CI, high CI
	tp.range=array(rep(NA,Nsp*2)) #range of plotTvals occupied
	dim(tp.range)=c(Nsp,2)
	ts.range=array(rep(NA,Nsp*2)) #range of species tvals
	dim(ts.range)=c(Nsp,2)
	bound=rep(0,Nsp*2) #boolean if species hits lower or upper boundary
	dim(bound)=c(Nsp,2)
	
	#loop through species and calculate niche breadth and intraspecific slopes
	# some of the if statements here apply to datasets where there is no data
	# for a species - this will normally only occur during bootstrapping
	for (i in minsp:maxsp) {
		#print(i)
		Ts = T[T$species==names(sp.N)[i],]
		if (nrow(Ts)>0) {
			tp.range[i,] = range(Ts$tp,na.rm=TRUE)
			ts.range[i,] = range(Ts$tval,na.rm=TRUE)
			if (tp.range[i,1]==plotT.r[1]) bound[i,1]=1
			if (tp.range[i,2]==plotT.r[2]) bound[i,2]=1
			sp.R[i] = diff(tp.range[i,])
			sp.Ntval[i]=length(Ts$tval[!is.na(Ts$tval)])
			
			if (calc.slopes) {
				if (use.spmeans) cases = 1:nrow(Ts) else 						cases = which(!is.na(Ts$tvalNA))
				ncases[i] = length(cases)
				if ((length(cases) > 1)&(sp.R[i]>0)) {
					#print(Ts)
					breg = summary(lm(Ts$tval[cases]~
						Ts$tp[cases],weights=Ts$abund[cases]))
					b[i,] = breg$coeff[,1]
					if (dim(breg$coeff)[1]>1) {
						b.se[i] = breg$coeff[2,2]
						b.p[i] = breg$coeff[2,4]
					}
					if (ncases[i] == 2) b.p[i] = 1
					#if (is.na(coeff[i,2])&!is.na(coeff[i,1])) 						#coeff[i,2]=0
					ureg = summary(lm(Ts$tval[cases]~Ts$tp[cases]))
					u[i,] = ureg$coeff[,1]
					if (dim(ureg$coeff)[1]>1) {
						u.se[i] = ureg$coeff[2,2]
						u.p[i] = ureg$coeff[2,4]
					}
					if (ncases[i] == 2) u.p[i] = 1
					#if(is.na(coeff[i,2])
						#&!is.na(coeff[i,1])) coeff[i,2]=0
					# Standardized major axis analysis
					# Remove comments here to run this analysis
					# Change method to 'MA' or 'OLS' for
					# alternative analyses
					#sm[i,] = 
					#	as.numeric(line.cis(Ts$tval[cases],					#	Ts$tp[cases],
					#	method='SMA')[2,])
				}
			}
		}
	}
		
	S <- data.frame(rownames(sp.N),sp.N,sp.abund,sp.T,
		sp.betaT,sp.alphaT,sp.R,ncases,b,b.se,b.p,u,u.se,u.p,
		tp.range,ts.range, bound)
	names(S) = c('species','Nplots','mean.abund','ts','betaT',
		'alphaT','Rs','slopeN','b.inter','bs','bs.se','bs.p',
		'u.inter','us','us.se','us.p',
		'min.tp','max.tp','min.tval','max.tval',
		'low.bound','hi.bound')
	return(S)
	#to obtain results of standardized major axis regression,
	#use alternative outputs below, and comment out those above
	#S <- data.frame(rownames(sp.N),sp.N,sp.abund,sp.T,
	#	sp.betaT,sp.alphaT,sp.R,ncases,b,b.se,b.p,u,u.se,u.p,
	#	sm,tp.range,ts.range, bound)
	#names(S) = c('species','Nplots','mean.abund','ts','betaT',
	#	'alphaT','Rs','slopeN','b.inter','bs','bs.se','bs.p',
	#	'u.inter','us','us.se','us.p',
	#	'sm.slope','sm.025','sm.975',
	#	'min.tp','max.tp','min.tval','max.tval',
	#	'low.bound','hi.bound')
	#return(S)
	
}
	
#Species jacknife
#this function removes each species in succession, calculates
#plot parameters, and then calculates parameters of that species
#against independently measured plot values; returns S matrix with
#same variables as above
#set compare = TRUE to view series of diagnostic plots comparing
#  values from S matrix based on normal vs. jackknife analysis
#comparing results of jackknife to full analysis.
calc.spec.attr.jackknife = function(D,T,S,weighted,
	calc.slopes=TRUE,use.spmeans=TRUE,compare=FALSE) {
	spp.list = unique(D$species)
	Nsp = length(spp.list)
	Sj = S
	for (i in 1:Nsp) {
		Tspp = T
		spp = spp.list[i]
		dj = D[-which(D$species==spp),]
		Tj = calc.means(dj,weighted)
		tpj = tapply(Tj$tp,Tj$plot,mean)
		Tspp$tp = 
			tpj[match(Tspp$plot,names(tpj))]
		Sspp = calc.spec.attr(Tspp,i,i,calc.slopes=calc.slopes,
			use.spmeans=use.spmeans)
		Sj[i,]=Sspp[i,]
		#plotTS(Sspp,Tspp)
	}
	#use compare=TRUE to examine S parameter values calculated
	#using all species vs. by species jackknife
	if (compare) {
		par(mfrow=c(2,1))
		par(ask=TRUE)
		for (j in 3:(ncol(S)-4)) {
			plot(S[,j],Sj[,j],main=names(S)[j])
			abline(0,1)
		}
		par(ask=FALSE)
	}
	return(Sj)
}

# Calculate plot attributes
# T has 7 columns: species, site, abund, tval, tp, ts, tvalNA
# S has 2 columns: see above
# This function assumes no missing values; NAs for traits should be filled in with species means.
#Returns data frame with 12 columns:
#plot	plot number of name
#Nsp	number of species
#tp	mean trait value
#Trange	range of trait values
#mnNiBr	mean niche breadth (S$Rs) of species in plot
#mnIntraSlope	mean value of S$bs for species in plot
#min.aT	minimum value of S$alphaT of species in plot
#$max.aT	maximum value of S$alphaT of species in plot
#tval.low	minimum trait value of species in plot
#tval.hi	maximum trait value of species in plot
#lowbound.sp	plot contains species with S$min.bound=1*
#hibound.sp	plot contains species with S$max.bound=1*


calc.plot.attr <- function(S,T) {
	plot.N <- tapply(T$tval, as.factor(T$plot), length)

	plot.T = tapply(T$tp,T$plot,median)

	Nplots <- length(plot.N)
	plotT.r <- range(plot.T)	
		
	bound.sp <- array(rep(0,Nplots*2))
	dim(bound.sp) <- c(Nplots,2)
	inplot.r <- array(rep(NA,Nplots*2))
	dim(inplot.r) <- c(Nplots,2)
	plotTdiv <- array(rep(NA,Nplots))
	plot.nibr <- array(rep(NA,Nplots))
	plot.Int <- array(rep(NA,Nplots))
	plot.dTr <- array(rep(NA,Nplots*2))
	dim(plot.dTr) <- c(Nplots,2)
	
	for (i in 1:Nplots) {
		#print(i)
		Tp = T[T[,2]==rownames(plot.N)[i],]
		sp.in.plot = S$species%in%Tp[,1]
		inplot.r[i,] = range(Tp$tval,na.rm=TRUE)
		plotTdiv[i] = diff(inplot.r[i,])
		plot.nibr[i] = mean(S$Rs[sp.in.plot],na.rm=TRUE)
		plot.Int[i] = mean(S$bs[sp.in.plot],na.rm=TRUE)
		plot.dTr[i,] = range(S$alphaT[sp.in.plot],na.rm=TRUE)
		if (sum(S$low.bound[sp.in.plot])>0) bound.sp[i,1]=1
		if (sum(S$hi.bound[sp.in.plot])>0) bound.sp[i,2]=1
		}
	P=data.frame(rownames(plot.N),plot.N,plot.T,
		plotTdiv,plot.nibr,plot.Int,plot.dTr,inplot.r,bound.sp)
	names(P) = c('plot','Nsp','tp','Trange',
		'mnNiBr','mnIntraSlope','min.aT','max.aT',
		'tval.low','tval.hi',
		'lowbound.sp','hibound.sp')
	return(P)
	}

#modified resampling function that returns x if length x = 1
resample <- function(x, size, replace=TRUE)
  if(length(x) <= 1) { if(!missing(size) && size == 0) x[FALSE] else rep(x,size)
  } else sample(x, size, replace=replace)

#boot.strata returns a bootstrapped sample of vector x,
#where bootstrapping is conducted independently within each level
#of variable strata. R's boot function allows for strata, but does not
#have the option of returning the actual bootstrapped sample
boot.strata = function(strata) {
	ns = length(strata)
	index = 1:ns
	strata = tapply(1:ns,as.numeric(strata))
	levels = length(unique(strata))
	for (i in 1:levels) {
		#print(i)
		rows = which(strata==i)
		index[rows] = resample(index[rows],	
			size=length(rows),replace=TRUE)
	}
	#print(index)
	return(index)
}

#bootstrapD conducts bootstrap of D matrix and calculates bootstrap values
#of betaT, alphaT, Rs, bs and us, outputting means, standard
#errors, and confidence intervals for each species. For slopes, the number of
#successful reps is also output, because for small sample
#sizes it is not possible to calculate a slope from all
#bootstrap samples. Returns a list. Element 1 records whether
#bootstrapping was stratified or not, and number of reps; Element 2 
#is a data frame with a row for each species, and columns for 
#statistics above. 
#
#If stratified = TRUE, then bootstrap is conducted independently
#within each species, so number of occurrences per species will
#be maintained, but number of species per plot will vary; if FALSE, 
#then all rows of T are resampled with replacement, and number of 
#occurrences per species will vary (sometimes = 0 for species with
#small sample sizes). In bootstrapping T, the plot mean trait
#values are resampled along with species trait values, rather than
#being recalculated from the bootstrap sample. The idea is that
#the original data gives an estimate of the distribution and 
#association between trait values and plot
#mean trait values of the kinds of plots in which the species
#lives, and we want confidence intervals on the resulting 
#parameters. If the original D matrix is sent to this function
#and line 10-11 are uncommented, line 12 commented out, then 
#you can bootstrap D and recalculate the resulting plot
#mean trait values as part of the resampling process.
bootstrapD = function(D,T,reps=1000,stratified=TRUE,
		showplot=FALSE,showreps=FALSE,
		calc.slopes=TRUE,use.spmeans=TRUE) {
	Tr = nrow(T)
	
	
	tval.bs = NULL
	beta.bs = NULL
	alpha.bs = NULL
	Rs.bs = NULL
	bs.bs = NULL
	us.bs = NULL
	Rab = NULL
	for (i in 1:reps) {
		if (showreps) print(i)
		if (stratified) br = 	
			boot.strata(as.numeric(T$species)) 			else br = sample(1:Tr,replace=TRUE)
		dbs = D[br,] #uncomment these two lines to bs D
		Tbs = calc.means(dbs,weighted)
		#Tbs = T[br,] #comment out this line to bs datain
		Sbs = calc.spec.attr(Tbs,
			calc.slopes=calc.slopes,use.spmeans=use.spmeans)
		tval.bs = rbind(tval.bs,Sbs$ts)
		beta.bs=rbind(beta.bs,Sbs$betaT)
		alpha.bs=rbind(alpha.bs,Sbs$alphaT)
		Rs.bs = rbind(Rs.bs,Sbs$Rs)
		bs.bs=rbind(bs.bs,Sbs$bs)
		us.bs=rbind(us.bs,Sbs$us)
		if (showplot) plotTS(Sbs,Tbs,weighted)
		Rab = rbind(Rab,cor(Sbs$betaT,Sbs$alphaT,use='pair'))
	}	
	bs.res = summary.BSstats(tval.bs,beta.bs,alpha.bs,Rs.bs,
		bs.bs,us.bs)
	bs.res = cbind(Sbs$species,bs.res)
	names(bs.res)[1]='species'
	return(list(c(as.character(stratified),reps),
		bs.res))
}

#Bootstrap two D matrices using same sampling for each
#to get distribution of bivariate trait, betaT 
#and alphaT correlations. Returns data frame
#with row for each bootstrap rep and columns for
#correlations of trait values, beta values, and
#alpha values.
bootstrap2D = function(D1,D2,reps=1000,stratified=TRUE,
		showplot=FALSE,weighted=TRUE,showreps=TRUE,
		calc.slopes=TRUE,use.spmeans=TRUE) {
	d1r = nrow(D1)
	d2r = nrow(D2)
	if (d1r!=d2r) {
		print('D matrices of different size')
		break
	}
	
	Rtt = NULL
	Rbb = NULL
	Raa = NULL
	for (i in 1:reps) {
		if (showreps) print(i)
		if (stratified) br = 	
			boot.strata(as.numeric(D1$species)) 			else br = sample(1:d1r,replace=TRUE)
		D1bs = D1[br,]
		D2bs = D2[br,]
		T1bs = calc.means(D1bs,weighted)
		T2bs = calc.means(D2bs,weighted)
		S1bs = calc.spec.attr(T1bs,calc.slopes=calc.slopes,
			use.spmeans=use.spmeans)
		S2bs = calc.spec.attr(T2bs,calc.slopes=calc.slopes,
			use.spmeans=use.spmeans)
		par(mfrow=c(1,3))
		if (showplot) {
			plot(S1bs$ts,S2bs$ts,cex=2,pch=19)
			plot(S1bs$betaT,S2bs$betaT,cex=2,pch=19)
			plot(S1bs$alphaT,S2bs$alphaT,cex=2,pch=19)
		}
		Rtt = rbind(Rtt,cor(S1bs$ts,S2bs$ts,use='pair'))
		Rbb = rbind(Rbb,cor(S1bs$betaT,S2bs$betaT,use='pair'))
		Raa = rbind(Raa,cor(S1bs$alphaT,S2bs$alphaT,use='pair'))
	}	
	return(data.frame(Rtt,Rbb,Raa))
}


#NULL MODEL of community assembly
#This null model randomizes elements of the input matrix, 
#breaking up observed associations among species, plots, trait 
#values, etc. The 'randomize' function may take as arguments one
# or more of the following:
#'plot', 'species', 'abund', 'traits', and 'abund+traits'.
#The corresponding columns will be randomized, and the data
#reanalyzed using the trait-gradient method.
#Ackerly and Cornwell paper uses: randomize plot: as a result, 
#species and their associated abundance and
#trait values are randomly assigned to plots; species mean trait
#values are maintained; expectation of beta is equal for all
#species, with wide range of alpha values
#
#nullD returns a list with three entries
#[[1]]: randomize option and number of reps
#[[2]]: dataframe with a row for each species and columns:
#			nullN: number of null reps + 1
#			tval: mean trait value under null
#			tval.se: standard error of null distribution of tval
#			tval.025: lower 2.5% CI of null distribution
#			tval.975: upper 97.5% CI of null distribution
#			betaT/se/025/975: beta value, se and 95% CI under null
#			alphaT/se/025/975: alpha value, se and 95% CI
#			Rs/se/025/975: niche breadth, se and 95% CI
#			bs.N: number of reps where slope could be calculated
#			bs/se/025/975: abundance weighted regression slope,
#				se and 95% CI
#			us.N: number of reps where unwtd slope was calculated
#			us/se/025/975: unweighted regression slope,
#				se and 95% CI
#[[3]]: dataframe with various statistics under null hypothesis
#			row 1 is observed data, subsequent rows are null reps
#			columns:
#			Rbeta-alpha: correlation of spp. beta and alpha values
#			beta.range: range of species beta values
#			alpha.range: range of species alpha values
#			Rs.range: range of niche breadths
#			tp.range: range of plot mean trait values


nullD = function(D,S,T,reps=999,
		randomize = c('plot'),showplot = FALSE,	
		calc.slopes=FALSE,use.spmeans=use.spmeans,showreps=TRUE) {
	dr = nrow(D)
	
	#Assign observed results to first rep of null
	ts.null = S$ts
	beta.null = S$betaT
	alpha.null = S$alphaT
	Rs.null = S$Rs
	bs.null = S$bs
	us.null = S$us
	stats = c(cor(S$betaT,S$alphaT,use='pair'),
			diff(range(S$betaT,na.rm=TRUE)),
			diff(range(S$alphaT,na.rm=TRUE)),
			diff(range(S$Rs,na.rm=TRUE)),
			diff(range(T$tp,na.rm=TRUE)))

	#null: assign traits at random to plots
	for (i in 1:reps) {
		if (showreps) print(i)
		dnull = D
		if ('plot' %in% randomize) 
			dnull$plot = D$plot[sample(1:dr,replace=FALSE)]
		if ('species' %in% randomize) 
			dnull$species = 	
				D$species[sample(1:dr,replace=FALSE)]
		if ('abund' %in% randomize) 
			dnull$abund = D$abund[sample(1:dr,replace=FALSE)]
		if ('traits' %in% randomize) 
			dnull[,4:5] = datain[sample(1:dr,replace=FALSE),4:5]
		if ('abund+traits' %in% randomize) 
			dnull[,3:5] = datain[sample(1:dr,replace=FALSE),3:5]
		Tnull = calc.means(dnull,weighted)
		Snull = calc.spec.attr(Tnull,
			calc.slopes=calc.slopes,use.spmeans=use.spmeans)
		ts.null = rbind(ts.null,Snull$ts)
		beta.null=rbind(beta.null,Snull$betaT)
		alpha.null=rbind(alpha.null,Snull$alphaT)
		Rs.null = rbind(Rs.null,Snull$Rs)
		bs.null=rbind(bs.null,Snull$bs)
		us.null=rbind(us.null,Snull$slopeNA)
		Rab = cor(Snull$betaT,Snull$alphaT,use='pair')
		tp.range = diff(range(Tnull$tp,na.rm=TRUE))
		beta.range = diff(range(Snull$betaT,na.rm=TRUE))
		alpha.range = diff(range(Snull$alphaT,na.rm=TRUE))
		Rs.range = diff(range(Snull$Rs,na.rm=TRUE))
		stats = rbind(stats,
			c(Rab,beta.range,alpha.range,Rs.range,tp.range))
		if (showplot) plotTS(Snull,Tnull,weighted)
	}	
	null.res = summary.BSstats(ts.null,beta.null,alpha.null,
		Rs.null,bs.null,us.null)
	names(null.res)[1]='nullN'
	stats = data.frame(stats)
	names(stats) = c('Rbeta-alpha','beta.range','alpha.range',
		'Rs.range','tp.range')
	return(list(c(randomize,reps),null.res,stats))
}

#null2D conducts null models for two traits to test the
#significance of correlations between trait parameters

null2D = function(D1,D2,S1,S2,T1,T2,reps=999,
		randomize = c('plot'),showplot = TRUE,showreps=TRUE
		,calc.slopes=FALSE) {
	d1r = nrow(D1)
	d2r = nrow(D2)
	if (d1r!=d2r) {
		print('T matrices of different size')
		break
	}
		
	#Assign observed results to first rep of null
	Rtt = cor(S1$ts,S2$ts,use='pair')
	Rbb = cor(S1$betaT,S2$betaT,use='pair')
	Raa = cor(S1$alphaT,S2$alphaT,use='pair')
	stats = c(Rtt,Rbb,Raa)

	#null: assign traits at random to plots
	for (i in 1:reps) {
		if (showreps) print(i)
		D1N = D1
		D2N = D2
		nsamp = sample(1:d1r,replace=FALSE)
		if ('plot' %in% randomize) {
			D1N$plot = D1$plot[nsamp]
			D2N$plot = D2$plot[nsamp]	
		}
		if ('species' %in% randomize) {
			D1N$species = D1$species[nsamp]
			D2N$species = D2$species[nsamp]	
		}
		if ('abund' %in% randomize) {
			D1N$abund = D1$abund[nsamp]
			D2N$abund = D2$abund[nsamp]	
		}
		if ('traits' %in% randomize) {
			D1N[,4:5]= D1[nsamp,4:5]
			D2N[,4:5]= D2[nsamp,4:5]
		}
		if ('abund+traits' %in% randomize) {
			D1N[,3:5]= D1[nsamp,3:5]
			D2N[,3:5]= D2[nsamp,3:5]
		}
		T1N = calc.means(D1N,weighted)
		T2N = calc.means(D2N,weighted)
		S1N = calc.spec.attr(T1N,calc.slopes=calc.slopes)
		S2N = calc.spec.attr(T2N,calc.slopes=calc.slopes)
		Rtt = cor(S1N$ts,S2N$ts,use='pair')
		Rbb = cor(S1N$betaT,S2N$betaT,use='pair')
		Raa = cor(S1N$alphaT,S2N$alphaT,use='pair')
		stats = rbind(stats,c(Rtt,Rbb,Raa))
		#if (showplot) plotTS(Snull,Tnull,weighted)
	}	
	stats = data.frame(stats)
	names(stats) = c('Rtt','Rbb','Raa')
	return(stats)
}

												
												
boxplot.one.null = function(N) {
	N = as.matrix(N)
	nr = nrow(N)
	nc = ncol(N)
	dim(N) = c(nr*nc,1)
	Nnam2 = rep(names(N),c(nr,nr,nr,nr,nr))
	N = data.frame(Nnam,as.numeric(N))
	obs = c(1,1+nr,1+2*nr,1+3*nr,1+4*nr)
	boxplot(N[,2]~N[,1])
	points(N[obs,1],N[obs,2],cex=2,pch=19,col='red')
	}		
	
	
#summarize statistics of bootstrap or null model samples
summary.BSstats = function(tval,beta,alpha,Rs,bs,us) {
	bs.res = NULL
	for (s in 1:ncol(beta))  		
	{
		bsN = length(which(!is.na(beta[,s])))
		tval.mn = mean(tval[,s],na.rm=TRUE)
		tval.se = sd(tval[,s],na.rm=TRUE)
		tval.ci = quantile(tval[,s],c(0.025,0.975),na.rm=TRUE)
		beta.mn = mean(beta[,s],na.rm=TRUE)
		beta.se = sd(beta[,s],na.rm=TRUE)
		beta.ci = quantile(beta[,s],c(0.025,0.975),na.rm=TRUE)
		alpha.mn = mean(alpha[,s],na.rm=TRUE)
		alpha.se = sd(alpha[,s],na.rm=TRUE)
		alpha.ci = quantile(alpha[,s],c(0.025,0.975),na.rm=TRUE)
		Rs.mn = mean(Rs[,s],na.rm=TRUE)
		Rs.se = sd(Rs[,s],na.rm=TRUE)
		Rs.ci = quantile(Rs[,s],c(0.025,0.975),na.rm=TRUE)
		bs.N = length(which(!is.na(bs[,s])))
		bs.mn = mean(bs[,s],na.rm=TRUE)
		bs.se = sd(bs[,s],na.rm=TRUE)
		bs.ci = quantile(bs[,s],c(0.025,0.975),na.rm=TRUE)
		us.N = length(which(!is.na(us[,s])))
		us.mn = mean(us[,s],na.rm=TRUE)
		us.se = sd(us[,s],na.rm=TRUE)
		us.ci = quantile(us[,s],
			c(0.025,0.975),na.rm=TRUE)
		bs.res = rbind(bs.res,c(bsN,tval.mn,tval.se,tval.ci,
			beta.mn,beta.se,beta.ci,						alpha.mn,alpha.se,alpha.ci,Rs.mn,Rs.se,Rs.ci,
			bs.N, bs.mn, bs.se,bs.ci,
			us.N, us.mn, us.se, us.ci))
	}
	bs.res=data.frame(bs.res)
	names(bs.res) = c('bsN',
		'tval','tval.se','tval.025','tval.975',
		'betaT','betaT.se','betaT.025','betaT.975',
		'alphaT','alphaT.se','alphaT.025','alphaT.975',
		'Rs','Rs.se','Rs.025','Rs.975',
		'bs.N','bs','bs.se','bs.025','bs.975',
		'us.N','us','us.se',
		'us.025','us.975')
	return(bs.res)
}


		
#extract list of species present in a plot
#plot can be identified by its mean trait value (will pickup
#all plots with identical value) or by plot name
extract.spec.from.plots=function(P,T,G=-99,pname='zzzzz') {
	if (G!=-99) plotname=P[P$tp==G,1] else plotname=pname
	Tp = T[T$plot==plotname,]
	return(Tp$species)
	}
	
#This plots local T vs. plot mean, and then add lines by species
#that have N > 3
plotTS <- function(S,T,spec.means=TRUE,slopes=TRUE,	linewd=1,minN=3,ylimits=NULL,linecol='red',
	pointcol='black',trait='trait value') {
		xlabel=paste('plot mean',trait)
		ylabel=paste('species',trait)
		plot(T$tp,T$tval,xlab=xlabel,ylab=ylabel,
			cex.lab=1.5,ylim=ylimits,pch=19,col=pointcol)
		Nsp = length(S$species)
		for (i in 1:Nsp) {
			if (S$Nplots[i] >= minN) {
				miny <- S$b.inter[i] + S$bs[i]*S$min.tp[i]
				maxy <- S$b.inter[i] + S$bs[i]*S$max.tp[i]
				if(slopes) lines(x=c(S$min.tp[i],
					S$max.tp[i]),
					y=c(miny,maxy),col=linecol,lwd=linewd)
				}
			}
		if (spec.means) points(S$betaT,S$ts,pch=19,cex=2)
		lmall = lm(T$tval~T$tp,weights=T$abund)
		abline(lmall,lty=2,lwd=2)
	}
	
plot.species <- function(S,T,snum,spec.means=TRUE,
	linewd=1,scol=c('red','blue','green','yellow'),
		ptype=c(15,16,17,8),mtype=c(22,21,24,8),
		pointcol='gray',trait='trait value') {
		xlabel=paste('plot mean',trait)
		ylabel=paste('species',trait)
		plot(T$tp,T$tval,xlab=xlabel,ylab=ylabel,
			cex.lab=1.5,pch=19,col=pointcol)
		Nsp = length(S$species)
		for (j in 1:length(snum)) {
			if (class(snum)=='numeric') i = snum[j]
			if (class(snum)=='character') 
				i = as.numeric(rownames(S)[snum[j]==S[,1]])
			if (S$Nplots[i] > 0) {
				points(T[T$species==S[i,'species'],'tp'],					T[T[,'species']==S[i,'species'],'tval'],
					pch=ptype[j],col=scol[j],cex=1.5)
				miny <- S$b.inter[i] + S$bs[i]*S$min.tp[i]
				maxy <- S$b.inter[i] + S$bs[i]*S$max.tp[i]
				lines(x=c(S$min.tp[i],S$max.tp[i]),
					y=c(miny,maxy),col=scol[j],lwd=linewd)
				if (spec.means) points(S$betaT[i],
					S$ts[i],pch=mtype[j],cex=2)
				}
			}
		abline(0,1,lty=2)
	}

write.tg.files = function(dname,D,T,S,Sna,P) {
	if(file.exists(dname)) {
		print('Directory exists - continue (1=Y)')
		reply=scan()
		if (reply==1) {
			setwd(dname)
			file.remove(list.files())
			setwd('..')
			file.remove(dname)
			} else break
		}
	dir.create(dname)
	setwd(dname)
	write.csv(D,file='D.csv',quote=FALSE,row.names=FALSE)
	write.csv(T,file='T.csv',quote=FALSE,row.names=FALSE)
	write.csv(S,file='S.csv',quote=FALSE,row.names=FALSE)
	write.csv(Sna,file='Sna.csv',quote=FALSE,row.names=FALSE)
	write.csv(P,file='P.csv',quote=FALSE,row.names=FALSE)
	setwd('..')
	}
	
read.tg.files = function(dname) {
	dexists=file.exists(dname)
	if (!dexists) {
		dir.create(dname)
		}
	setwd(dname)
	D=read.csv('D.csv')
	T=read.csv('T.csv')
	S=read.csv('S.csv')
	Sna=read.csv('Sna.csv')
	P=read.csv('P.csv')
	setwd('..')
	return(list(D,T,S,Sna,P))
	}