FDvar=function(abun,trait){#a single-column matrix of abundances; a matrix of traits
  abunNoNA=as.numeric(abun[!is.na(trait)])
  abunNoNA=abs(abunNoNA)
  probNoNA=abunNoNA/sum(abunNoNA)
  trait=trait[!is.na(trait)]
  lnxbar=sum(probNoNA*log(trait-min(trait)+1))#eqn2 in Mason 2003 JVS; forced traits to all positive values
  V.k=vector()
  for (k in 1:length(trait)){
    Vk=probNoNA[k]*((abs(log(trait[k]-min(trait)+1)-lnxbar))^2)
    V.k=c(V.k,Vk)
  }
  V=sum(V.k)
  return((2/pi)*atan(5*V))
}