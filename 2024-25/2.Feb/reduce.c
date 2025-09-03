
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  int myrank, size, sendval[3], minval[3], maxval[3], sumval[3];
  MPI_Status status;

  struct {
   float val;
   int rank;
  } buf, recvbuf;

  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  sendval[0] = myrank; // initialization
  sendval[1] = myrank+1; // initialization
  sendval[2] = myrank+2; // initialization

  // reduce count = 1 
  MPI_Reduce(sendval, maxval, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); // find max
								       
  // reduce count = 3 
  MPI_Reduce(sendval, minval, 3, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD); // find min of sendvals
  MPI_Reduce(sendval, sumval, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // find sum of sendvals

  if (!myrank) printf ("reduced %d\n%d %d %d\n%d %d %d\n", maxval[0], minval[0], minval[1], minval[2], sumval[0], sumval[1], sumval[2]); //random check

  // different root
  
  float flnum= 2.67*(myrank+1), val, val_, val_s, val_p;

  MPI_Reduce(&flnum, &val_s, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&flnum, &val_s, 1, MPI_FLOAT, MPI_SUM, 5, MPI_COMM_WORLD); 

  MPI_Reduce(&flnum, &val_p, 1, MPI_FLOAT, MPI_PROD, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&flnum, &val_p, 1, MPI_FLOAT, MPI_PROD, 5, MPI_COMM_WORLD); 

  if (myrank== 0 || myrank == 5) 
     printf ("%d: %f %f\n", myrank, val_s, val_p);

  MPI_Finalize();
  return 0;

}

