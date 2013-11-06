### Code for Q3 : Analysis of the breast cancer data ####

setwd("/Users/rohosenb/Documents/UCDavis/STA-250/HW1/mcmcgibbs") ###setting directory

####Including the libraries I am gonna need
library(MASS)
library(coda)

####### Function calculating the logarithm of the posterior distribution ###############
logpost <- function(beta,m,y,X,beta.0,Sigma.0.inv){
  ans = sum(y*(X%*%beta)) - 0.5*t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0) - sum(m*log(1 + exp(X%*%beta)))
  return(ans)
}

################ 1st Method: Implementing function bayes.logreg --- It will do an MH of dimension "size" #########################
bayes.logreg <- function(m,size,y,X,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=T){
  beta.in <- rep(0,size) ####Initial beta
  post.sample <- NULL ####Sample of beta from its posterior distribution to be produced
  beta <- beta.in
  sig <- diag(size)   #### variance of proposal density
  ac <- 0
  for(i in 1:(niter+burnin)){
    if(i <= burnin){
      temp <- mvrnorm(1,beta,sig)
      u <- runif(1)
      if(logpost(temp,m,y,X,beta.0,Sigma.0.inv) - logpost(beta,m,y,X,beta.0,Sigma.0.inv) > log(u))
      {
        beta <- temp
        ac <- ac + 1
      }
      
      #####printing acceptance rate at every retune iterations
      if(i%%retune==0){
        print(paste(c("acceptance rate at step ", i,"is =", ac*100/retune,"%"),sep=" "))
        ac <- 0
      }
      post.sample <- rbind(post.sample,beta)
    }
    else{
      temp <- mvrnorm(1,beta,sig)
      u <- runif(1)
      if(logpost(temp,m,y,X,beta.0,Sigma.0.inv) - logpost(beta,m,y,X,beta.0,Sigma.0.inv) > log(u)){
        beta <- temp
      }
      post.sample <- rbind(post.sample,beta)
    }
    
    ####### Printing the values of beta at every print.every-th iteration
    if(i%%print.every==0) print(paste(c("beta value at iteration",i,"=",beta),sep=" "))
  }
  return(post.sample)
}

########################## 2nd Method:  Using MH within Gibbs ############################################################

################## MH within Gibbs ################################
bayes.metro <- function(m,size,y,X,beta_init,beta.0,Sigma.0.inv){
  sig <- diag(size)
  temp <- mvrnorm(1,beta_init,sig)
  u <- runif(1)
  if(logpost(temp,m,y,X,beta.0,Sigma.0.inv) - logpost(beta_init,m,y,X,beta.0,Sigma.0.inv) > log(u)) beta <- temp
  else beta <- beta_init
  return(beta)
}

######################## Gibbs sampling #############################
bayes.gibbs <- function(m,y,X,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE){
  beta <- rep(sample(1:50,1),11)
  post.sample <- NULL
  for(i in 1:(burnin+niter)){
    for(j in 1:11) beta[j] <- bayes.metro(m=1,size=1,y,X=as.matrix(X[,j]),as.matrix(beta[j]),as.matrix(beta.0[j]),as.matrix(Sigma.0.inv[j,j]))	
    post.sample <- rbind(post.sample,beta)
  }
  return(post.sample)
}


################################ Applying the methods on the breast cancer dataset ##################################
###### Set up spcification ###################
cancer = read.table("cancer.txt",header=T)
y = (cancer$diagnosis=="M")        #### Response 
X = as.matrix(cancer[,1:10])
Xs = cbind(1,scale(X))             #### Design Matrix
beta.0 <- matrix(0,11,1)
Sigma.0.inv <- diag(11)/1000


########################### Without Gibbs --- 11-dimensional MH on the Cancer data : Implementing the 1st method ################################
res_MH = bayes.logreg(m=1,size=11,y,X=Xs,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100)

######Plots
pdf("MH plots.pdf")
par(mfrow=c(4,3))
for(i in 1:11) plot(res_MH[,i],type="l",main=paste("beta_",i-1,sep=""))
dev.off()

##### Calculating effective sample-size
MC_MH = as.mcmc(res_MH[1001:11000,])    ### Converting my posterior sample matrix into an MCMC object
effectiveSize(MC_MH)


########## With Gibbs ---- Using the conditional distribution of beta_j given the rest -- MH within Gibbs : Implementing the 1st method###########

res_Gibbs = bayes.gibbs(m=1,y,X=Xs,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE)

######Plots
pdf("Gibbs plots.pdf")
par(mfrow=c(4,3))
for(i in 1:11) plot(res_Gibbs[1001:11000,i],type="l",main=paste("beta_",i-1,sep=""))
dev.off()

##### Calculating effective sample-size
MC_Gibbs = as.mcmc(res_Gibbs[1001:11000,])
effectiveSize(MC_Gibbs)

### I compared the two methods, MH within Gibbs seems to work lot better than the 11-dimensional MH. 
### So I used the sample I obtained from the 2nd method, for the rest of the analysis

######## Lag 1 autocorrelation #########
auto=NULL
for(i in 1:11) auto[i] = cor(res_Gibbs[1001:10999,i],res_Gibbs[1002:11000,i])

######### 95% Credible Interval ###########
int = matrix(0,11,2)
for(i in 1:11) int[i,] = matrix(quantile(res_Gibbs[1001:11000,i],probs=c(0.025,0.975)),1,2)

###### Gelman Rubin #########
for(i in 1:5){
  ch = bayes.gibbs(m=1,y,X=Xs,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE)
  chain[[i]] = as.mcmc(ch)
}
chain = mcmc.list(chain)
gelman.diag(chain)
  
########### Diagnostic ###############
###Parameter generation###
beta_true = mvrnorm(100,beta.0,solve(Sigma.0.inv))
chap = Xs%*%t(beta_true)
acid1 = rep(0,11)
acid2 = 0
for(i in 1:100){
  ###Data generation###
  
  data = matrix(rbinom(569,1,c(exp(chap[,i])/(1+exp(chap[,i])))))
  
  ###Simulate from posterior###
  check = bayes.gibbs(m=1,data,X=Xs,beta.0,Sigma.0.inv,niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE)
  
  ###Credible interval###
  int_check = matrix(0,11,2)
  for(j in 1:11) int_check[j,] = matrix(quantile(check[1001:11000,j],probs=c(0.025,0.975)),1,2)
  
  ### Check for inclusion ###
  for(j in 1:11) acid1[j] = acid1[j] + ((beta_true[i,j] >= int_check[j,1])&(beta_true[i,j] <= int_check[j,2]))
  include = 1
  for(j in 1:11) include = include * ((beta_true[i,j] >= int_check[j,1])&(beta_true[i,j] <= int_check[j,2]))
  acid2 = acid2 + include
}
write.table(c(acid1,acid2),file="coverage.txt")`

################# Done -- happy! :) #########################