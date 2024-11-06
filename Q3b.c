#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#ifdef _OPENMP
#include<omp.h>
#endif

double q(int i, int j, double delta){
  double x = -1.0 + i*delta;
  double y = -1.0 + j*delta;
  return 2.0*(2.0 - (x*x) - (y*y));
}


int main(int argc, char* argv[])
{
  double delta;
  int thread_count;
  if (argc == 3)
    {
      delta = atof(argv[1]);
      thread_count = strtol(argv[2],NULL,10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }
  
  int N = (int) 2.0/delta + 1;
  
  double phi[N][N];
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      phi[i][j] = 0.0;
    }
  }
  
  double exact_phi[N][N];
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      double x = -1 + i*delta;
      double y = -1 + j*delta;
      exact_phi[i][j] = (x*x -1)*(y*y -1);
    }
  }

  float error_percent(double exact_phi[N][N], double phi[N][N], int N){
    float error_num = 0;
    float error_den = 0;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < N; j++){
	error_num += (phi[i][j] - exact_phi[i][j])*(phi[i][j] - exact_phi[i][j]);
	error_den += exact_phi[i][j]*exact_phi[i][j];
      }
    }
    return (float) error_num / error_den;
  }
  
  int iter = 0;
  float error = 1.0;
  int nd = 2*(N-1) - 1;
  int i_begin, i_end;
  int i, j, l;
  
  double start_time, end_time;
  start_time = omp_get_wtime();
  
  /*#pragma omp parallel default(none) num_threads(thread_count) private(i,j) shared(phi,l, i_begin, i_end, delta, error, N, nd, exact_phi,iter)*/
  while(error > 0.01){
  for(l = 0; l < nd; l++){
    if(l < N-1){
      i_begin = 0;
      i_end = l+1;
    }
    else{
      i_begin = l-N+2;
      i_end = N-1;
    }
  
    /*#pragma omp for*/
 #pragma omp parallel for default(none) num_threads(thread_count) private(i,j) shared(phi,l, i_begin, i_end, delta, N)
    for(i = i_begin; i < i_end; i++){
      j = l-i+1;
      if(i == 0 || j == 0 || i == N-1 || j == N-1){
	phi[i][j] = 0.0;
      }
      else{
	phi[i][j] = ( phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] + (delta*delta*q(i,j,delta)) )*0.25;
      }
    }
  }
  error = error_percent(exact_phi, phi, N);
  printf("%d \n",iter);
  iter++;
  }
 end_time = omp_get_wtime();
 
 printf("time taken: %lf \n",end_time - start_time);
 printf("number of iterations: %d \n",iter);
 
 j = (int) 1.5/delta ;
  for(int i = 0; i < N; i++){
    printf("%f %f %f \n", -1.0 + i*delta, exact_phi[i][j], phi[i][j]);
     }

}
