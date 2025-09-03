#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[])
{
  
  int myrank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  int rows = atoi(argv[1]);
  int cols = atoi(argv[2]);

  int color = myrank%cols;
  int newrank, newsize;
  MPI_Comm newcomm;

  //split into sub-communicators based on color (odd-even in this example)
  MPI_Comm_split (MPI_COMM_WORLD, color, myrank /*preserving original rank order*/, &newcomm);

  MPI_Comm_rank( newcomm, &newrank );
  MPI_Comm_size( newcomm, &newsize );

  //MPI_Gather ();

  MPI_Comm_free(&newcomm);

  printf("Rank %d: new rank %d of %d\n", myrank, newrank, newsize);

  MPI_Finalize();
  return 0;

}
