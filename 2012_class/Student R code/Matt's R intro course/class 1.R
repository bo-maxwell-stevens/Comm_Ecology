# SECTION 2.2 Inputing data into R (c, scan)
# c creates a numerical data vector (1 dimension) series of numbers or letters
c(1,2,3)
# TO SAVE THIS SERIES of numbers YOU NEED TO MAKE AN OBJECT
x=c(1,2,3)
abc=c("a","b","c")
# scan allows you to input values sequentially
scan()
1
2
3
x

#need a blank line to end
# can copy paste a column of values from excel and use scan() to read them in
# this can be useful to store data in simple test eg 1 sample-t
# can use vectors to do other things than importing data, eg creating loops (repetitions)
#but usually you will have more than just one set of data so...

# SECTION 2.3: Data Types (numeric, character, matrix array, data.frame, list)
# use MATRIX command
matrix()
#matrix takes 1-D data and turns it into a rowsxcolumns format
x=c(1,2,3,4,5,6,7,8,9)
matrix(x)
matrix(x,nrow=3,ncol=3)
matrix(nrow=3,ncol=3)
m=matrix(x,nrow=3,ncol=3)
m

#another way to make a matrix using the ARRAY command. can make n-dimensional matrices!
array(data=x, dim=c(3,3))
#dim=c() tell how many rows and col (r,c) there are going to be
array(data=x,dim=c(2,3,2))
#array(data=x,y) then array(scan(,y)) if you have already specified an x and y. then list all the x's using scan
# instead of listing numbers, eg c(1,2,3,4,5), you can say (1:5) meaning numbers 1-5
m=array(x,c(3,3))
m
#the bracets in matrices give where you are, eg row and col. useful for data manipulation
m[2,2]
#this will give the value at row2 col2
m[2,c(2,3)]
m[2,c(2:3)]
#this gives values of row2 for col 2 & 3
#need the c bc if you dont have it, it will read m[2,2,3] which reads like 3 dimensions (2nd row, 2nd col, 3rd dim).
colnames(m)=c("sp.a","sp.b","sp.c")
rownames(m)=c("obs.d","obs.e","obs.f")
m
#this is for NUMERICAL data only. to put in characters in the matrix eg male and felmale, need to use another data object called the DATA FRAME

df<-data.frame(m)
m
#to change 1st col into male/female
df[,1]<-c("m","f","m")
#to run this, you can just highlight the df and puch command return
df[1,1]<-"f"
#this overwrites the previous designation of [1,1] from male to female from the other df

#if you want to add on another column
df=data.frame(m)
sex=c("m","f","m")
cbind(df,sex)

#take one object eg df and add something to it. the order is important. cbind is for column, rbind is for row
rbind(df,sex)

log.sp.a<-log(df[,1])
df<-cbind(df,log.sp.a)
df


# How to know what type of data you have: is. fct
#is.matrix returns a logical (true or false). is this a matrix?
is.matrix(m)
is.matrix(df)

# class fct: what the data type is
class(df)
#this will pull up data.frame

our.list<-list(m,df)
#list is conglomeration of mult matrices into one big list. to go into the list use double bracets
our.list[[2]][2,2]
#this gets the data in the 2nd object in the list from row2 col2

# if you want to detatch a package use: detach(package:______)