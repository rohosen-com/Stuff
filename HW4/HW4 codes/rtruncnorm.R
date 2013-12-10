library(RCUDA)
library(truncnorm)

#setting the values
N = 10000L
mu = rep(2, N)
sigma = rep(1, N)
a = rep(0, N)
b = rep(1.5, N)
rng_a = 89
rng_b = 56
rng_c = 1
maxtries = 5000
verbose = TRUE

if (verbose){
    cat("Setting cuGetContext(TRUE)...\n")
}
cuGetContext(TRUE)

if (verbose){
    cat("Loading module...\n")
}
m <- loadModule("rtruncnorm.ptx")
if (verbose){
    cat("done. Extracting kernel...\n")
}
rtruncnorm_kernel <- m$rtruncnorm_kernel

# Note: kernel looks like this:
# __global__ void rtruncnorm_kernel(float *x, int n, 
#                          float *mu, float *sigma, 
#                          float *a, float *b,
#                          int rng_a, int rng_b, int rng_c,
#                          int maxtries);

if (verbose){
    cat("done. Setting up miscellaneous stuff...\n")
}


### generate random numbers using function rtruncnorm
x = rtruncnorm(N, a, b, mean = mu, sd = sigma)

threads_per_block <- 512L
block_dims <- c(threads_per_block, 1L, 1L)
grid_d1 <- as.integer(ceiling(N/threads_per_block))
grid_dims <- c(grid_d1, 1L, 1L)

if (verbose){
    cat("Grid size:\n")
    print(grid_dims)
    cat("Block size:\n")
    print(block_dims)
}

nthreads <- prod(grid_dims)*prod(block_dims) 
if (verbose){
    cat("Total number of threads to launch = ",nthreads,"\n")
}
if (nthreads < N){
    stop("Grid is not large enough...!")
}

cudaDeviceSynchronize()

if (verbose){
    cat("Running CUDA kernel...\n")
}
cu_time <- system.time({
    if (verbose){
        cat("Copying truncated random N(mu, sigma)'s to device...\n")
    }
    cu_copy_to_time <- system.time({mem <- copyToDevice(x)})
    # .cuda(kernel, args, gridDim, blockDim)
    cu_kernel_time <- system.time({.cuda(rtruncnorm_kernel, mem, N, mu, sigma, a, b, rng_a, rng_b, rng_c, maxtries, 
                                         gridDim=grid_dims, blockDim=block_dims, outputs=NULL)})
    if (verbose){
         cat("Copying result back from device...\n")
    }
    cu_copy_back_time <- system.time({cu_ret <- copyFromDevice(obj=mem,nels=mem@nels,type="float")})
    # Equivalently:
    #cu_ret <- mem[]
})


cu_kernel2_time <- system.time({y <- .cuda(rtruncnorm_kernel, "x"=x, N, mu, sigma, a, b, rng_a, rng_b, rng_c, maxtries, 
                                           gridDim=grid_dims, blockDim=block_dims, outputs="x")})

cat("")
r_time <- system.time({
    r_ret <- rtruncnorm(N, a, b, mean = mu, sd = sigma)
})
if (verbose){
    cat("done. Finished profile run! :)\n")
}
# Not the best comparison but a rough real-world comparison:
if (verbose){
    cat("CUDA time:\n")
    print(cu_time)
    cat("R time:\n")
    print(r_time)
    cat("Breakdown of CUDA time:\n")
    cat("Copy to device:\n")
    print(cu_copy_to_time)
    cat("Kernel:\n")
    print(cu_kernel_time)
    cat("Copy from device:\n")
    print(cu_copy_back_time)
    cat("Smart kernel version timing:\n")
    print(cu_kernel2_time)
}

# Differences due to floating point vs. double...
tdiff <- sum(abs(cu_ret - r_ret))
if (verbose){
    cat("Diff = ",tdiff,"\n")
    cat("Differences in first few values...\n")
    print(abs(diff(head(cu_ret)-head(r_ret))))
    cat("Differences in last few values...\n")
    print(abs(diff(tail(cu_ret)-tail(r_ret))))
}

# Experiment with clean up:
if (verbose){
    cat("Device memory usage prior to cleanup:\n")
    print(cuMemGetInfo())
    rm(mem)
    gc()
    cat("Device memory usage prior after cleanup:\n")
    print(cuMemGetInfo())
}




