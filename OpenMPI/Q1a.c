#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#define PI 3.141593

double u_x(double x){
  return sin(4*PI*x); 
}

double initial_conditions(int i, double delx, double(*u_x)(double)){
  if (delx*i <= 0.5)
    {
      return u_x(delx*i);
    }
  else{
    return 0.0;
  }
  
}

double upwind_scheme(int i, double delx, double delt, double u_n[], double c){
  double gamma;
  gamma = c*(delt/delx);
  
  return u_n[i]*(1 - gamma) + gamma*u_n[i-1];
  
}

double QUICK_scheme(int i, double delx, double delt, double u_n[], double c){
  double gamma;
  gamma = c*(delt/delx);
  
  return u_n[i]*(1.0 - (3.0/8.0)*gamma) + gamma*((7.0/8.0)*u_n[i-1] - (1.0/8)*u_n[i-2] - (3.0/8)*u_n[i+1]);
  
}


int main(int argc, char* argv[])
{
  double L = 2.0;
  double delx = 0.002;
  double delt = 0.0001;
  double c = 1.0;
  double t;
  if (argc == 2)
    {
      t = atof(argv[1]);
      
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }
  
  int N = (int)L/delx;
  
  double output[N];

  clock_t t_start;
  clock_t timetaken;
  t_start = clock();
  
  output[0] = 0;
  output[N-1] = 0;
  
  for (int i = 1; i <= N-2  ; i++){
    output[i] = initial_conditions(i, delx, u_x);
  }

  
  if (t != 0.0){
    for (int j = 0; j < t/delt; j++){
      
      for (int i = 1; i <= N-2; i++){

	/*if(i == 1){*/
         output[i] = upwind_scheme(i, delx, delt, output, c);
	 
	 /*}
	
      else{
	  output[i] = QUICK_scheme(i, delx, delt, output, c);
	  
	  }*/
       }
    }
}
  
 timetaken = clock() - t_start;
 
 printf("time taken: %.7f \n",(float)timetaken/CLOCKS_PER_SEC);
 printf("\n");
 
 for(int i = 0; i < N; i ++){
   printf(" %.3f \n",output[i]);
   }

 printf("\n");

    }
