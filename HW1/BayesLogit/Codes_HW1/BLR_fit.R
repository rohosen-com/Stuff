##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(mvtnorm)
library(coda)

########################################################################################
########################################################################################
## Handle batch job arguments:



# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################



################ Fit the Bayesian Logistic Regression model to the dataset#########################

library(MASS)

data = read.table(file=paste("data/blr_data_", sim_num, ".csv", sep=""),header=TRUE, sep=",")

####### Function calculating the logarithm of the posterior distribution ###############
logpost <- function(beta,m,y,X,beta.0,Sigma.0.inv)
{
  ans = sum(y*(X%*%beta)) - 0.5*t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0) - sum(m*log(1 + exp(X%*%beta)))
  return(ans)
}

########### Function calculating all the percentiles of a sample #####################
percent=c(1:99)/100
percentile <- function(sample)  ###sample is a vector of the sample of which percentiles you want to calculate
{
  percentiles = quantile(sample,percent,names=FALSE)
  return(percentiles)
}

####################### Implementing function bayes.logreg#########################################
bayes.logreg <- function(m, y, X, beta.0, Sigma.0.inv, niter, burnin, print.every,retune,verbose)
{
  
  beta.in <- c(10,10)  ####Initial beta
  post.sample <- NULL  ####Sample of beta from its posterior distribution to be produced
  beta <- beta.in
  sig <- diag(2)  #### variance of proposal density
  ac <- 0
  
  ########### I consider a proposal density of N(0,diag(1,1))
  
  ###################Burnin part ###################################
  for(i in 1:burnin)
  {
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
    
    ####### Printing the values of beta at every print.every-th iteration
    if(i%%print.every==0) print(paste(c("beta value at iteration",i,"=",beta),sep=" "))
    post.sample <- rbind(post.sample,beta)
  }
  
  ######################### Actual iteration part######################
  for(i in (burnin+1):(burnin+niter))
  {
    temp <- mvrnorm(1,beta,sig)
    u <- runif(1)
    if(logpost(temp,m,y,X,beta.0,Sigma.0.inv) - logpost(beta,m,y,X,beta.0,Sigma.0.inv) > log(u)){
      beta <- temp
    }
    post.sample <- rbind(post.sample,beta)
    if(i%%print.every==0) print(paste(c("beta value at iteration",i,"=",beta),sep=" "))
  }
  
  
  
  
  ########## Calculating the percentiles of the samples drawn from the posterior distribution using MH #########
  
  percentiles.beta=matrix(0,nrow=99,ncol=2)
  percentiles.beta[,1]=percentile(post.sample[(burnin+1):(burnin+niter),1])
  percentiles.beta[,2]=percentile(post.sample[(burnin+1):(burnin+niter),2])
  
  write.table(percentiles.beta, file=paste("results/blr_res_", sim_num, ".csv",sep=""), row.names = FALSE, col.names = FALSE, sep=",")

}

######################### Using the function bayes.logreg on my dataset ##########################
##Set up spcification
m <- as.matrix(data$n)
y <- as.matrix(data$y)
X <- cbind(data$X1,data$X2)
beta.0 <- matrix(0,2,1)
Sigma.0.inv <- diag(2)

bayes.logreg(m, y, X, beta.0=c(0,1), Sigma.0.inv, niter=10000, burnin=1000, print.every=1000,retune=100,verbose=TRUE)


cat("done. :)\n")
