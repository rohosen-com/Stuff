####################### RCUDA code for Probit MCMC ##################################
library(RCUDA)
library(MASS)
library(truncnorm)
setwd("/Users/rohosenb/Documents/UCDavis/STA-250/HW4")

### Read data
data = read.table("data_01.txt", header = TRUE)
data = as.matrix(data)

#### Set  values
n = as.integer(nrow(data))
y = data[,1]                    # response y
X = data[, -1]                  # design matrix X

### Prior values
beta_0 = rep(0,ncol(X))
Sigma_0_inv = matrix(0, length(beta_0), length(beta_0))

### No. of iterations
niter = 2000
burnin = 500

#### Calculating sum(x_i*x_i) for a matrix X
summies <- function(X,n){
  sum=0
  for(i in 1:n){
    sum = sum + X[i,]%*%t(X[i,])
  }
   return(sum) 
}
XtX = summies(X,n)

######### Figuring out the truncation range for Z ###########
##### Truncation range (a,b)
  a = rep(0, n)
  b = rep(0, n)
  a[which(y == 0)] = -Inf # y= 0 means (-inf, 0)
  b[which(y == 1)] = Inf  # y=1 means (0, inf)

##################### Initial values of Z ########################
Z_int = rnorm(n) ### generate Z from N(0,1) as that's distribution of Z under beta_0 = 0
### Check if value of Z matches with y or change
y[which(y == 0)] = -1
Zy =Z_int*y
Z_int[which(Zy < 0)] = (-1)*Z_int[which(Zy < 0)]


##### Generating sample from the Conditonal posterior distribution of beta given Z #######
beta_sample <- function(z, X, XtX, beta_0, sigma_0_inv ){
  
  ### Calculating sum of z_i*x_i
  zx = z*X
  sum_zx = colSums(zx)
  
  ### Setting the mu and sigma of conditional distribution of beta
  sigma_beta = solve(sigma_0_inv + XtX)
  mu_beta = sigma_beta%*%(sigma_0_inv%*%beta_0 + sum_zx)
  
  #### Sampling beta from multivariate N(mu_beta, sigma_beta)
  beta_sample = mvrnorm(1,mu_beta,sigma_beta)
  beta_sample = as.vector(beta_sample)
  return(beta_sample)
}

############################ The Probit MCMC function for CPU #############################################
probit_mcmc_cpu = function(
  n,           # length of y, Z
  Z_int,       # initial value of Z 
  X,           # (n x p) design matrix
  XtX,         # enumerated XX' 
  a,           # lower bound of truncation 
  b,           # upper bound of truncation
  beta_0,      # (p x 1) prior mean
  Sigma_0_inv, # (p x p) prior precision 
  niter,       # number of post burnin iterations
  burnin      # number of burnin iterations
)
{
  Z = Z_int
  sd = rep(1,n)
  beta.sample = NULL
  
  for(t in 1:(burnin+niter)){
    beta = beta_sample(Z, X, XtX, beta_0, Sigma_0_inv)
    beta.sample = rbind(beta.sample, beta)   ### the columns of beta.sample gives the different betas from different iterations
    
    mean = X%*%beta
    Z = rtruncnorm(1, a, b, mean, sd)
  }
  betas.sample = as.matrix(beta.sample[(burnin+1): (burnin+niter), ])
  beta_sample_mean = colMeans(betas.sample)
  write.table(beta_sample_mean, file=paste("posterior_mean_beta_CPU_", n, ".txt",sep=""), 
              row.names = FALSE, col.names = FALSE, sep=" ")
}

#### Run the function and measure the system time
CPU.time = system.time(probit_mcmc_cpu(n, Z_int, X, XtX, a, b, beta_0, Sigma_0_inv, niter, burnin))

#### print the system time
cat("CPU time:\n")
print(CPU.time)



####################### The Probit MCMC function for GPU ########################################################

########## Setting up the block and grid size 
threads_per_block <- 512L
block_dims <- c(threads_per_block, 1L, 1L)
grid_d1 <- as.integer(ceiling(n/threads_per_block))
grid_dims <- c(grid_d1, 1L, 1L)

##### The probit MCMC function

probit_mcmc_gpu = function(
  n,           # length of y, Z
  Z_int,       # initial value of Z 
  X,           # (n x p) design matrix
  XtX,         # enumerated XX' 
  a,           # lower bound of truncation 
  b,           # upper bound of truncation                
  beta_0,      # (p x 1) prior mean
  Sigma_0_inv, # (p x p) prior precision 
  niter,       # number of post burnin iterations
  burnin,      # number of burnin iterations
  block_dims,  # block dimensions
  grid_dims    # grid_dimensions
  )

{
    rng_a = 89
    rng_b = 56
    maxtries = 1000
    beta.sample.GPU <- NULL
    Z.GPU = Z_int
    sd = rep(1,n)
    
    cuGetContext(TRUE)
    m <- loadModule("rtruncnorm.ptx")
    rtruncnorm_kernel <- m$rtruncnorm_kernel
    
    nthreads <- prod(grid_dims)*prod(block_dims) 
    if (nthreads < n){
      stop("Grid is not large enough...!")
    }
    
    #### Staring the iterations to sample beta and Z's from the coditionals
    
    for(t in 1:(burnin+niter)){
      
      #### Set rng_c to iteration no.
      rng_c = t
      
      ### Draw sample of beta using initial value of Z
      beta = beta_sample(Z.GPU, X, XtX, beta_0, Sigma_0_inv)
      beta.sample.GPU = rbind(beta.sample.GPU, beta)   ### the columns of beta.sample gives the different betas from different iterations
      
      ### Draw sample of Z using beta(t)
      # mean of Z using current beta(t)
      mean = X%*%beta
      
      cudaDeviceSynchronize()
      
      .cuda(rtruncnorm_kernel, "Z.GPU"=x, n, mean, sd, a, b, rng_a, rng_b, rng_c, maxtries, 
            gridDim=grid_dims, blockDim=block_dims, outputs="Z.GPU")})
      
      
    }

    ##### Write the betas in an output file 
    betas.sample.GPU = as.matrix(beta.sample.GPU[(burnin+1): (burnin+niter), ])
    beta_sample_mean.GPU = colMeans(betas.sample.GPU)
    write.table(beta_sample_mean.GPU, file=paste("posterior_mean_beta_GPU_", n, ".txt",sep=""), 
            row.names = FALSE, col.names = FALSE, sep=" ")
  }

#### Run the function and measure the system time
GPU.time = system.time(probit_mcmc_gpu(n, Z_int, X, XtX, a, b, beta_0, Sigma_0_inv, niter, burnin, block_dims, grid_dims))

#### print the system time
cat("GPU time:\n")
print(GPU.time)






