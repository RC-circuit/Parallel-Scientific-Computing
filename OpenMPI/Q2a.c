#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include <stdbool.h>
#define PI 3.141593

double q(int i, int j, double delta){
  double x = -1.0 + i*delta;
  double y = -1.0 + j*delta;
  return (x*x) + (y*y);
}

double error(int N, double phit_1[][N], double phit[][N])
  {
    double norm = 0.0;
    for(int j = 0; j < N; j++){
      for(int i = 0; i < N; i++){
        norm += pow((phit[i][j] - phit_1[i][j]),2);
    }
  }
  return sqrt(norm);
  }

int main(int argc, char* argv[])
{
  double delta = 0.1;
  int N = (int) 2.0/delta + 1;
  
  double phit[N][N];
  double phit_1[N][N];
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      if(i == 0){
	phit[i][j] = sin(2*PI*(-1 + j*delta));
	phit_1[i][j] = sin(2*PI*(-1 + j*delta));
      }
      else{
        phit[i][j] = 0;
        phit_1[i][j] = 0;
      }
    }
  }
  
  
  int iteration = 0;

  clock_t start_time, end_time;
  double time_taken;
  bool condition = true;
  start_time = clock();
  while(condition){
  for(int i = 1; i <= N-1; i++){
    for(int j = 1; j < N-1; j++){
      if(i == N-1){
	phit_1[i][j] = (4.0*phit_1[i-1][j] - phit_1[i-2][j])/3.0;
      }
      else{
      phit_1[i][j] = ( phit[i+1][j] + phit[i-1][j] + phit[i][j+1] + phit[i][j-1] + (delta*delta*q(i,j,delta)) )*0.25;
      }
    }
  }
  condition = (error(N, phit_1, phit) > 0.0001);
  iteration ++;
  
  for(int j = 0; j < N; j++){
    for(int i = 0; i < N; i++){
	phit[i][j] = phit_1[i][j];
    }
  }
  
  }
  
  end_time = clock();
  printf("time taken %lf \n", (double) (end_time - start_time)/CLOCKS_PER_SEC);
  printf("number of iterations: %d \n",iteration);

  
  /*for(int j = 0; j < N; j++){
    printf("%f \n",phit[j][N/2]);
    }*/
  
}
