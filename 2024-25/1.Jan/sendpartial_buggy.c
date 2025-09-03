#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  MPI_Status status;
  double sTime, eTime, time;

  MPI_Init(&argc, &argv);

  // count is the number of elements desired to be sent
  // count_ is the total number of elements in the array
  int count = atoi (argv[1]), count_ = atoi (argv[2]);
  int buf[count];
  int buffer[count_];
  int size, recvcount=0, myrank; 

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // initialize data
  for (int i=0; i<count_; i++)
   buf[i] = myrank+i;

  if (myrank == 0)
   MPI_Send (buf, count, MPI_INT, 1, 1, MPI_COMM_WORLD);	   
  else
   if (myrank == 1) {
    MPI_Recv (&buf[count], count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	   
  }

  printf ("%d %d %d\n", myrank, buf[9], buf[10]);

  MPI_Finalize();
  return 0;

}

