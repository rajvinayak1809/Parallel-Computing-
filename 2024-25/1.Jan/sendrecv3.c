// Timing codes

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size; 
  MPI_Status status;
  double sTime, eTime, time;

  MPI_Init(&argc, &argv);

  int count = atoi (argv[1]);
  int buf[count];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // initialize data
  for (int i=0; i<count; i++)
   buf[i] = myrank+i;

  if (myrank == 0) {
    MPI_Send (buf, count, MPI_INT, 1, 1, MPI_COMM_WORLD);	   
    MPI_Recv (buf, count, MPI_INT, 1, myrank, MPI_COMM_WORLD, &status);	   
  }
  else
   if (myrank == 1) {
    MPI_Send (buf, count, MPI_INT, 0, 0, MPI_COMM_WORLD);	   
    MPI_Recv (buf, count, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);	   
   }

  printf ("%d %d %d\n", myrank, count, status.MPI_SOURCE);

  MPI_Finalize();
  return 0;

}

