// Timing codes

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, buf[100000];
  int count = atoi(argv[1]);

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  // initialize data
  for (int i=0; i<count; i++)
      buf[i] = myrank + i*i;

  // has to be called by all processes
  double sTime = MPI_Wtime();
  MPI_Bcast(buf, count, MPI_INT, 1, MPI_COMM_WORLD);
  double eTime = MPI_Wtime();

  // simple check
  printf ("%d: %lf %d %d %d\n", myrank, eTime - sTime, buf[0], buf[1], buf[2]);

  MPI_Finalize();
  return 0;

}

