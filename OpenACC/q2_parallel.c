#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE		float
#define SMALLVALUE	0.001

void initialise(int N, int ngangs, TYPE mat[][N])
{
  #pragma acc parallel loop collapse(2) present(mat[0:N][0:N]) num_gangs(ngangs)
  for (int ii = 0; ii < N; ++ii)
    for (int jj = 0; jj < N; ++jj)
      {	mat[ii][jj] = (ii + jj) / (float)N / N;}

  #pragma acc parallel loop present(mat[0:N][0:N]) num_gangs(ngangs)
  for (int ii = 0; ii < N; ++ii)
    mat[ii][ii] = 1.0;
}
			
			
void printMat(int N, int ngangs, TYPE a[][N])
{
  #pragma acc serial present(a[0:N][0:N]) num_gangs(ngangs)
  {
  for (int ii = 0; ii < N; ++ii)
    {
      for (int jj = 0; jj < N; ++jj)
	printf("%.2f ", a[ii][jj]);
      printf("\n");
    }
  }
}

void cholesky(int N, int ngangs, TYPE a[][N])
{
  #pragma acc parallel loop present(a[0:N][0:N]) num_gangs(ngangs)
  for (int ii = 0; ii < N; ++ii) {
    #pragma acc loop gang
    for (int jj = 0; jj < ii; ++jj) {
      #pragma acc loop worker
      for (int kk = 0; kk < jj; ++kk){
	    a[ii][jj] += -a[ii][kk] * a[jj][kk];
      }
      a[ii][jj] /= (a[jj][jj] > SMALLVALUE ? a[jj][jj] : 1);
      
    }
  }
    #pragma acc parallel loop present(a[0:N][0:N]) num_gangs(ngangs)
    for (int ii = 0; ii < N; ++ii) {
        #pragma acc loop gang
        for (int kk = 0; kk < ii; ++kk){
            a[ii][ii] += -a[ii][kk] * a[ii][kk];
            }
        a[ii][ii] = sqrt(a[ii][ii]);
    }
  
}

int main(int argc, char* argv[])
{
    
    int N;
    int ngangs;

    if (argc == 3)
    {
      N = atoi(argv[1]);
      ngangs = atoi(argv[2]);     
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }

  TYPE a[N][N];

  struct timespec start, end;
  

  #pragma acc data create (a[0:N][0:N])
  { 
  
  clock_gettime(CLOCK_MONOTONIC, &start);

  initialise(N, ngangs, a);
  cholesky(N, ngangs, a);
  //printMat(N, ngangs, a);

  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / 1e9;
  printf("Elapsed time: %f seconds\n", elapsed);

  }

  
  return 0;
}
