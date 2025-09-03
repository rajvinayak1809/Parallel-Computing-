#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, P, n, count;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &P);

  int myArraySize = atoi(argv[1]);
  double data[myArraySize], recvbuf[myArraySize];
  double time_, total;
  MPI_Status status;

   for (int i=0; i<myArraySize; i++) {
      data[i] = pow (2, i) * 20.01;
   }

  time_ = MPI_Wtime();

  if (P == 1 && P%2 != 0) return 0;

  if (myrank % 2 == 0 && myrank < P-1) // P-1 is anyway odd (assume)
  {
      // Send/recv right neighbor from even ranks
      MPI_Send (data, myArraySize, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
      MPI_Recv (recvbuf, myArraySize, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &status);
  }

  else if (myrank % 2 != 0 && myrank > 0) // 0 is anyway excluded here
  {
      // Send/recv left neighbor
      MPI_Recv (recvbuf, myArraySize, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &status);
      MPI_Send (data, myArraySize, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
  }

  if (myrank % 2 != 0 && myrank < P-1) // P-1 is anyway odd (assume)
  {
      // Send/recv right neighbor from odd ranks
      MPI_Send (data, myArraySize, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
      MPI_Recv (recvbuf, myArraySize, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &status);
  }

  else if (myrank % 2 == 0 && myrank > 0) // 0 is anyway excluded here
  {
      // Send/recv left neighbor
      MPI_Recv (recvbuf, myArraySize, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &status);
      MPI_Send (data, myArraySize, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
  }

  time_ = MPI_Wtime() - time_;

  // get the maximum time (slowest rank determines the overall time)
  MPI_Reduce(&time_, &total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) printf ("%lf\n", total);

  MPI_Finalize();
  return 0;

}

