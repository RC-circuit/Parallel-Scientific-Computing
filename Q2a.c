#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#ifdef _OPENMP
#include<omp.h>
#endif

float f_x(float x){
  return sin(5.0*x);
}

float exact(float x){
  return 5.0*cos(5.0*x);
}

int main(int argc, char* argv[])
{
  int  N;
  if (argc == 2)
    {
      N = strtol(argv[1],NULL,10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }
  
  float B[N+1];
  for(int i = 0; i < N+1; i++){
    if(i == 0 || i == N){
      B[i] = 1;
    }
    else{
      B[i] = 4;
    }
  }

  float A[N];
  for(int i = 0; i < N; i++){
    if(i == N-1){
      A[i] = 2;
    }
    else{
      A[i] = 1;
    }
  }

  float C[N];
  for(int i = 0; i < N; i++){
    if(i == 0){
      C[i] = 2;
    }
    else{
      C[i] = 1;
    }
  }

  float Y[N+1];
  float h = (float)3.0/N;
  for(int i = 0; i < N+1; i++){
    if(i == 0){
      Y[i] = ( (-2.5)*f_x(h*i) + 2.0*f_x(h*(i+1)) + 0.5*f_x(h*(i+2)) )*(1./h);
    }
    else if(i == N){
      Y[i] = (2.5*f_x(h*i) - 2.0*f_x(h*(i-1)) - 0.5*f_x(h*(i-2)) )*(1./h);
    }
    else{
      Y[i] = (f_x(h*(i+1)) - f_x(h*(i-1)))*(3./h);
    }
  }

  float U[N+1];
  float L[N];
  U[0] = B[0];
  for(int i = 0; i < N; i++){
    L[i] = A[i]/U[i];
    U[i+1] = B[i+1] - L[i]*C[i];
  }

float Z[N+1];
Z[0] = Y[0];
for(int i = 1; i < N+1 ; i++){
  Z[i] = Y[i] - L[i-1]*Z[i-1];
 }

 float X[N+1];
 X[N] = Z[N]/U[N];
 for(int i = N-1; i >= 0; i--){
   X[i] = (Z[i] - C[i]*X[i+1])*(1.0/U[i]);
 }

 for(int j = 0; j < N+1; j ++){
   printf("%f  %f  %f\n", h*j, X[j], exact(h*j));
  }
}
  
