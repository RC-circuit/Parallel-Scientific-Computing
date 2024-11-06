// OpenMP program to calculate first derivative using recursive double method
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

float func(float x)	{
	return (sin(5.0 * x));
}

float dfunc(float x)	{
	return(5.0 * cos(5.0 * x));
}

int main(int argc, char* argv[])	{
	
	int num_thread = 1;
	
	if(argc == 2)	{
		num_thread = strtol(argv[1], NULL, 10);
	}
	else	{
		printf("Command line argument of number of threads is required.\n");
		return 1;
	}
	
	int N = 1000;			// Number of divisions
	double Lx = 3.0;		// Length of domain
	
	// Allocate vectors
	float *xgrid = (float *)malloc((N+2) * sizeof(float));
	float *f = (float *)malloc((N+2) * sizeof(float));
	float *x = (float *)malloc((N+2) * sizeof(float));
	float *df_exact = (float *)malloc((N+2) * sizeof(float));
	float *a = (float *)malloc((N+2) * sizeof(float));
	float *b = (float *)malloc((N+2) * sizeof(float));
	float *c = (float *)malloc((N+2) * sizeof(float));
	float *d = (float *)malloc((N+2) * sizeof(float));
	float *alpha = (float *)malloc((N+2) * sizeof(float));
	float *beta = (float *)malloc((N+2) * sizeof(float));
	
	double h, t1, t2;
	int i, k, k_steps;
	FILE *fptr;

	h = Lx / N;	// grid spacing
	k_steps = ceil(log(N+1.0)/log(2.0));
	
	// Initialize vectors
	for(i = 1; i < (N+2); i++)	{
		xgrid[i] = (i-1) * h;
		f[i] = func(xgrid[i]);
		df_exact[i] = dfunc(xgrid[i]);	// Analytical solution
		x[i] = 0.0;
	}
	
	// Fill sub, diagonal, super and RHS vectors
	a[1] = 0.0;
	b[1] = 1.0;
	c[1] = 2.0;
	d[1] = (1.0 / h) * (-2.5 * f[1] + 2.0 * f[2] + 0.5 * f[3]);
	for(i = 2; i < (N+1); i++)	{
		a[i] = 1.0;
		b[i] = 4.0;
		c[i] = 1.0;
		d[i] = (3.0 / h) * (f[i+1] - f[i-1]);
	}
	a[N+1] = 2.0;
	b[N+1] = 1.0;
	c[N+1] = 0.0;
	d[N+1] = (1.0 / h) * (2.5 * f[N+1] - 2.0 * f[N] - 0.5 * f[N-1]);
	
	t1 = omp_get_wtime();
	#pragma omp parallel num_threads(num_thread) default(shared) private(i)
	{
	// Elimination phase
	for(k = 1; k <= k_steps; k++)	{
		#pragma omp for
		for(i = 1; i < (N+2); i++)	{
			if (i >= ((int)pow(2, k-1) + 1))	{
				alpha[i] = -a[i] / b[i - (int)pow(2, k-1)];	// Compute alpha[i]
			}
			else        {
				alpha[i] = 0.0;
			}
			
			if (i <= (N + 1 - (int)pow(2, k-1)))	{
				beta[i] = -c[i] / b[i + (int)pow(2, k-1)];	// Compute beta[i]
			}
			else        {
				beta[i] = 0.0;
			}
			
			if (i >= ((int)pow(2, k) + 1))	{
				a[i] = alpha[i] * a[i - (int)pow(2, k-1)];	// Compute a[i]
			}
			else        {
				a[i] = 0.0;
			}
			
			if (i <= (N + 1 - (int)pow(2, k)))	{
				c[i] = beta[i] * c[i + (int)pow(2, k-1)];	// Compute c[i]
			}
			else        {
				c[i] = 0.0;
			}
			
			b[i] = alpha[i] * c[i - (int)pow(2, k-1)] + b[i] + beta[i] * a[i + (int)pow(2, k-1)];	// Compute b[i]
			
			d[i] = alpha[i] * d[i - (int)pow(2, k-1)] + d[i] + beta[i] * d[i + (int)pow(2, k-1)];	// Compute d[i]
		}
	}
	}
	
	#pragma omp parallel num_threads(num_thread) default(shared) private(i)
	{
	#pragma omp for
	// Solution phase
	for(i = 1; i < (N+2); i++)	{
		x[i] = d[i] / b[i];
	}
	}
	t2 = omp_get_wtime();
	
	printf("Time taken for convergence = %lf\n", t2-t1);	
	
	// Deallocate memory
	free(xgrid);
	free(f);
	free(x);
	free(df_exact);
	free(a);
	free(b);
	free(c);
	free(d);
	free(alpha);
	free(beta);
		
	return 0;
}
