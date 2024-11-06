#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>

double f_x(double x){
  return sin(5.0*x);
}

double exact(double x){
  return 5.0*cos(5.0*x);
}

int main(int argc, char* argv[])
{
  int  N;
  int ngangs;
  
  if (argc == 3)
    {
      N = strtol(argv[1],NULL,10);
      ngangs = strtol(argv[2],NULL,10);
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }
  
  double A[N+1][N+1];
  double Y[N+1];
  double X[N+1];
  double h = (double)3.0/N;

  #pragma acc data create(A[0:N+1][0:N+1]) create(Y[0:N+1]) copyout(X[0:N+1]) copy(h)
  {

  #pragma acc parallel loop collapse(2) present(A[0:N+1][0:N+1])
  for(int i = 0; i < N+1; i++){
    for(int j = 0; j < N+1; j++){
        A[i][j] = 0.0;
    }
  }

  #pragma acc parallel loop present(A[0:N+1][0:N+1])
  for(int i = 0; i < N+1; i++){
    for(int j = 0; j < N+1; j++){
       if(i == j){
        if(i > 0 && i < N){
            A[i][i-1] = 1.0;
            A[i][i] = 4.0;
            A[i][i+1] = 1.0;
           }
        else{
         A[i][i] = 1.0;
           }
        }

        else if((i == 0) && (j == 1)){
        A[i][j] = 2.0;
            }
        else if((i== N) && (j == N-1)){
        A[i][j] = 2.0;
            }    
        }
    }

  #pragma acc parallel loop present(Y[0:N+1])
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
  
  for(int k = 0; k < N; k++){
    #pragma acc parallel loop present(A[0:N+1][0:N+1])
    for(int i = k+1; i < N+1; i++){
        A[i][k] = A[i][k]/A[k][k];
    }
    #pragma acc parallel loop collapse(2) present(A[0:N+1][0:N+1])
    for(int i = k+1; i < N+1; i++){
        for(int j = k+1; j < N+1; j++){
            A[i][j] -= A[i][k]*A[k][j];
        }
    }
  }


  
for(int k = 0; k < N ; k++){
    #pragma acc parallel loop num_gangs(ngangs) present(Y[0:N+1])
    for(int i = k+1; i < N+1; i++){
         Y[i] = Y[i] - A[i][k]*Y[k];
    }
 }


 double sum;
 
 #pragma acc serial present(X[0:N+1]) create(sum)
 {
  X[N] = Y[N]/A[N][N];
 for(int k = N-1; k >= 0; k--){
    sum = 0.0;
    for(int i = k+1; i < N+1; i++){
        sum += A[k][i]*X[i];
        }
    X[k] = (Y[k] - sum)/A[k][k];
    }    
 }

  }

//  for(int j = 0; j < N+1; j ++){
//    printf("%f  %f  %f\n", h*j, X[j], exact(h*j));
//   }

}
  