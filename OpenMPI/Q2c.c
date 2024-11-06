#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include <stdbool.h>
#include<mpi.h>
#define PI 3.141593

double q(int i, int j, double delta, int ln, int my_id){
  double x = -1.0 + (j*delta);
  double y = -1.0 + ((my_id*ln + i)*delta);
  return (x*x) + (y*y);
}

double error(int N, double phit_1[][N], double phit[][N], int ln)
  {
    double norm = 0.0;
    for(int i = 1; i < ln + 1; i++){
      for(int j = 0; j < N; j++){
        norm += pow((phit[i][j] - phit_1[i][j]),2);
    }
  }
  return norm;
  }


int main(int argc, char* argv[])
{

  int my_id, size, tag = 100;
  double L = 2.0;
  double delta;
  int iter = 0;
  int k;
  
  if (argc == 2) //if k is required to be taken as input change argc = 3
    {
      delta = atof(argv[1]);
      //k = atoi(argv[2]);
      
      
    }
  else
    {
      printf("\n A command line argument other than name of the executable is required...exiting the program..\n");
      return 1;
    }
  int N =  (int) L/delta + 1;

  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Barrier(MPI_COMM_WORLD);

  double start_time = MPI_Wtime();
  int ln = (int) (N / size);
  
  double phit[N];
  phit[N-1] = 0;
  double local_phit[ln+2][N];
  double local_phit_1[ln+2][N];
  double local_error;
  double global_error;
  bool condition = true;


  for(int i = 0; i < ln + 2; i++){
    for(int j = 0; j < N; j++){
        if(i == 0 || i == ln + 1){
            local_phit[i][j] = 0;
            local_phit_1[i][j] = 0;
        }
        else if(j == 0){
            local_phit[i][j] = sin(2*PI*(-1.0 + ((my_id*ln + i-1)*delta)));
            local_phit_1[i][j] = sin(2*PI*(-1.0 + ((my_id*ln + i-1)*delta)));
            }
        
        else{
            local_phit[i][j] = 0;
            local_phit_1[i][j] = 0;
        } 
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  while (condition){
    ////////////////send receive for even//////////////////////
  
  for (int i = 1; i < ln+1; i++){
    if ((my_id == 0 && i == 1)) {
        continue;
      }
    else{
        for (int j = 1; j < N-1 ; j++){
            if((i+j) % 2 == 0){
        local_phit_1[i][j] = (local_phit[i-1][j] + local_phit[i+1][j] + local_phit[i][j-1] + local_phit[i][j+1] + (pow(delta,2)*q(i, j, delta, ln, my_id)))*0.25;
        }
      }
      local_phit_1[i][N-1] = (4.0*local_phit_1[i][N-2] - local_phit_1[i][N-3])/3.0 ;
    }
  }

   if (my_id != size -1)
    {
      MPI_Send(&local_phit_1[ln][0], N, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD);
    }
  if (my_id != 0){
    MPI_Recv(&local_phit_1[0][0], N, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (my_id != 0){
    MPI_Send(&local_phit_1[1][0], N, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD);
  }
  if(my_id != size - 1){
    MPI_Recv(&local_phit_1[ln+1][0], N, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

////////////////send receive for odd//////////////////////
for (int i = 1; i < ln+1 ; i++){
      if ((my_id == 0 && i == 1)) {
        continue;
      }
    else{
        for (int j = 1; j < N-1 ; j++){
            if((i+j) % 2 == 1){
        local_phit_1[i][j] = (local_phit_1[i-1][j] + local_phit_1[i+1][j] + local_phit_1[i][j-1] + local_phit_1[i][j+1] + (pow(delta,2)*q(i, j, delta, ln, my_id)))*0.25;
        }
      }
      local_phit_1[i][N-1] = (4.0*local_phit_1[i][N-2] - local_phit_1[i][N-3])/3.0 ;
    }
}

if (my_id != size -1)
    {
      MPI_Send(&local_phit_1[ln][0], N, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD);
    }
  if (my_id != 0){
    MPI_Recv(&local_phit_1[0][0], N, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (my_id != 0){
    MPI_Send(&local_phit_1[1][0], N, MPI_DOUBLE, my_id - 1, tag, MPI_COMM_WORLD);
  }
  if(my_id != size - 1){
    MPI_Recv(&local_phit_1[ln+1][0], N, MPI_DOUBLE, my_id + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
 
 local_error = error(N, local_phit_1, local_phit, ln);

 MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 global_error = sqrt(global_error);
 condition = (global_error > 0.0001);
 iter++;
 
 MPI_Barrier(MPI_COMM_WORLD);

 for(int i = 0; i < ln+2; i++){
    for(int j = 0; j < N; j++){
        local_phit[i][j] = local_phit_1[i][j];
       }
    }
}

MPI_Barrier(MPI_COMM_WORLD);

double end_time = MPI_Wtime();

if(my_id == 0){
    printf("number of iterations: %d \n", iter);
    printf("time taken %f in sec\n",end_time - start_time);
  }

/*MPI_Barrier(MPI_COMM_WORLD);

if(my_id == size/2){
    for(int j = 0; j < N; j ++){
    printf("%f \n",local_phit[1][j]);    
       }
    }

MPI_Barrier(MPI_COMM_WORLD);

if(my_id == k){
for(int i = 1; i < ln+1; i++){
    printf("%f \n",local_phit[i][N/2]);
    }
}*/

  MPI_Finalize();

  return 0;
}