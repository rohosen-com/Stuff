n = 1000
k = 200
AR1gen = function(j){

#Q6.part(a): generating the AR(1)series for n
  e=rnorm(n)
  y=NULL 
y[1]=e[1]            #since y0=0, y1 is nothing but e1

#I am assuming roh_t is a const and its value is roh
#So the AR(1) I am simulating from is: y(t)=roh*y(t-1)+e(t)
roh=0.9
for (i in 2:n){
  y[i]=roh*y[i-1]+e[i]
}
  return(y)
}

data = sapply(1:k,AR1gen)

########################################################################
par(mfrow=c(2,2))
#taking mean of each of k realizations of the AR(1) process at each time-point
meanAR=rowMeans(data)
plot(meanAR,type="b",lwd=2,col="blue",pch=18, main="Mean plot at each time point",ylim=c(-1,1))

#taking Variance of each of k realizations of the AR(1) process at each time point
varAR=apply(data,1,var)
plot(varAR,type="b",lwd=2,col="red",pch=18, main="Variance plot at each time point",ylim=c(0,9))

#########################################################################
#taking mean of each of n realizations of the AR(1) process
meanAR1=colMeans(data)
plot(meanAR1,type="b",lwd=2,col="purple",pch=18, main="Mean plot of each AR(1) process",ylim=c(-1,1))

#taking Variance of each of n realizations of the AR(1) process
varAR1=apply(data,2,var)
plot(varAR1,type="b",lwd=2,col="orange",pch=18, main="Variance plot of each AR(1) process",ylim=c(0,9))


