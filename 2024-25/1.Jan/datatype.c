// Timing codes

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, recvcount=0; 
  MPI_Status status;
  double sTime, eTime, time;

  MPI_Init(&argc, &argv);

  int count = atoi (argv[1]), count_ = 10*count;
  int buf[count];
  char buffer[count_];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // initialize data
  for (int i=0; i<count; i++)
   buf[i] = myrank+i;

  if (myrank == 0)
   MPI_Send (buf, count, MPI_INT, 1, 1, MPI_COMM_WORLD);	   
  else
   if (myrank == 1) {
    MPI_Recv (buffer, count_, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);	   
    MPI_Get_count (&status, MPI_CHAR, &recvcount);
  }

  printf ("%d %d %d\n", myrank, count, recvcount);

  MPI_Finalize();
  return 0;

}

