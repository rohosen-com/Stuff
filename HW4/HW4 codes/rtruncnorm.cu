
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{

__global__ void rtruncnorm_kernel(float *x, int n, float *mu, float *sigma, float *a, float *b, int rng_a, int rng_b, int rng_c, int maxtries)
{

	/* variables needed */

	int accepted = 0;
	int righttrunc;
	int tries = 1;
	float u_exp, alpha, z, phi, u;

	/* Usual block/thread indexing... */
	int myblock = blockIdx.x + blockIdx.y * gridDim.x;
	int blocksize = blockDim.x * blockDim.y * blockDim.z;
	int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
	int idx = myblock * blocksize + subthread;

	/* Check if idx is less than n */
	if(idx > n){
	printf("No need to take samples anymore");
	return;
	}

	/* Initialize the random no. generator */
	curandState rng;
	curand_init(rng_a + idx*rng_b, rng_c, 0, &rng);
       
        /* Try the rejection sampling */

	while(!accepted && tries < maxtries){
		x[idx] = mu[idx] + sigma[idx]*curand_normal(&rng);
		if( x[idx] >= a[idx] && x[idx] <= b[idx]){
			accepted = 1;
			}
		tries = tries + 1;
		}
	
	/* If not successful use approximation from Robert(2009) */

	while(!accepted){

		/* Check whether left truncation or right truncation */

		if(isfinite(a[idx])) { 
			righttrunc = 1; 
			}
		else {
			righttrunc = 0;
			}

		/* If left truncation make it to a right one */
		if(!righttrunc){
			a[idx] = -b[idx];
			b[idx] = CUDART_NAN_F;
			}

		/********************************************************************/
		/* Step1: Generating a random number from exponential(alpha, a[idx]) */

		u_exp = curand_uniform(&rng);		
		alpha = (a[idx] + sqrt(a[idx]*a[idx] + 4))/2;
		z = a[idx] - (log(u_exp)/alpha);

		/* Step2: Computing phi */

		if(a[idx] < alpha){
			phi = exp(-(alpha - z)*(alpha - z)/2);
			}
		else {
			phi = exp(-(alpha - z)*(alpha - z)/2)*exp(-(a[idx] - alpha)*(a[idx] - alpha)/2);
			}

		/* Step3: Generating a random number from U[0,1] */
			
		u = curand_uniform(&rng);

		/* FinalStep: accepting or rejecting the random number */

		if(u < phi){ 
			accepted = 1;
			x[idx] = sigma[idx]*z + mu[idx];
			if(!righttrunc) x[idx] = sigma[idx]*(-z) + mu[idx];
			}

	}	
	return;
}

} // END extern "C"

























































