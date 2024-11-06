#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
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
  
  return u_n[i]*(1.0 - (3.0/8)*gamma) + gamma*((7.0/8)*u_n[i-1] - (1.0/8)*u_n[i-2] - (3.0/8)*u_n[i+1]);
  
}


int main(int argc, char* argv[])
{

  int my_id, size, tag = 100;
  double L = 2.0;
  double delx = 0.002;
  double delt = 0.0001;
  double c = 1.0;
  int N = (int)L/delx + 1;
  double output[N];
  double local_sum;
  
  
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
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  int ln = (int)((float)N/size);
  double local_output[ln + 1];

  output[0] = 0;
  output[N-1] = 0;
  
  if (my_id ==0){
  for (int i = 1; i <= N-2  ; i++){
    output[i] = initial_conditions(i, delx, u_x);
    }
  }

  MPI_Scatter(output, ln, MPI_DOUBLE, local_output + 1, ln, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  for (int j  = 0; j < t/delt; j++){
  if (my_id == 0)
    {
      MPI_Send(local_output + ln, 1, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD);
    }
  else if (my_id == size -1){
    MPI_Recv(local_output, 1, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_id % 2 == 0){
    MPI_Send(local_output + ln, 1, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD);
    MPI_Recv(local_output, 1, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else{
    MPI_Recv(local_output, 1, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(local_output +ln, 1, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD);
  }

  for (int i = 1; i < ln+1; i++){
    if ((my_id == 0 && i==1) || (my_id == size-1 && i==ln)) {
        continue;
      }

      else{
	local_output[i] = upwind_scheme(i, delx, delt, local_output, c);
      }
   }
 }

  double output_final[N];

  MPI_Gather(local_output + 1, ln, MPI_DOUBLE, output_final, ln, MPI_DOUBLE, 0 , MPI_COMM_WORLD);

  if (my_id == 0){
    for(int i = 0; i < N; i ++){
   printf(" %.3f \n",output_final[i]);
   }
  }
  MPI_Finalize();

  return 0;
  
}
  
