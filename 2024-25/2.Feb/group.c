#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]) 
{
  int numranks, rank, len;

  // initialize MPI
  MPI_Init (&argc, &argv);

  // get number of processes
  MPI_Comm_size (MPI_COMM_WORLD, &numranks);

  // get my rank
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  // get the world group
  MPI_Group g_group;
  MPI_Comm_group (MPI_COMM_WORLD, &g_group);

  // create new group ranks array for even ranks
  int ranks[numranks/2], i, j=-1;

  if (rank%2 == 0)
  {
   for (i=0; i<numranks; i+=2)
      ranks[++j] = i;
  }
  else  
  for (i=1; i<numranks; i+=2)
      ranks[++j] = i;

  int ranks_size = j+1; // generic way

  // create a new group
  MPI_Group new_group;
  MPI_Group_incl (g_group, ranks_size, ranks, &new_group);

  MPI_Comm new_comm;
  MPI_Comm_create_group (MPI_COMM_WORLD, new_group, 100, &new_comm);

  // must be called by all processes
  //MPI_Comm_create (MPI_COMM_WORLD, new_group, &new_comm);

  // size of new comm
  int new_size, new_rank;
  MPI_Comm_size (new_comm, &new_size);
  MPI_Comm_rank (new_comm, &new_rank);
  printf ("Old rank %d, new rank %d\n", rank, new_rank);

  // done with MPI
  MPI_Finalize();
}

