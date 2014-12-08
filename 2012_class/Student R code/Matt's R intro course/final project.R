#logistic growth equation:
# x(t+1)= x(t) + r*X(t)
x=1
r=2^x
z=x+x*r
y=0
for (i in 1:100) {y[i]=x[i]+x[i]*r}
y