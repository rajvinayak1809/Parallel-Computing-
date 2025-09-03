// Multiple senders 1 receiver
//
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( int argc, char *argv[])
{

  int BUFSIZE = atoi(argv[1]);
  int arr[BUFSIZE], recvarr[BUFSIZE];
  int myrank, size;
  MPI_Status status;
  double start_time;
  char *buf;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  for (int i=0;i<BUFSIZE; i++)
    arr[i] = myrank+i; 

  // buffer space
  int bufsize = BUFSIZE*5;
  buf = (char *) malloc(bufsize);

  // Only one buffer can be attached to a process at a time
  MPI_Buffer_attach(buf, bufsize);
  //MPI_Buffer_attach(buf, MPI_BUFFER_AUTOMATIC); //v4 onwards

  start_time = MPI_Wtime ();
  for (int j=0; j<1; j++) {
  if (myrank == 0){ 
     MPI_Bsend(arr, BUFSIZE/2, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD);
     MPI_Bsend(arr+BUFSIZE/2, BUFSIZE/2, MPI_INT, myrank+2, myrank+2, MPI_COMM_WORLD);
     //MPI_Bsend(arr+BUFSIZE/2, BUFSIZE/2, MPI_INT, myrank+2, myrank+2, MPI_COMM_WORLD);
  MPI_Buffer_detach(&buf, &bufsize);
  }
  else 
     MPI_Recv(recvarr, BUFSIZE, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
     //MPI_Recv(recvarr, BUFSIZE/2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  if (myrank == 0) printf ("Rank %d: Send time= %lf\n", myrank, MPI_Wtime() - start_time);
  else printf ("%d\n", recvarr[0]);
  }


  MPI_Finalize();
  return 0;
}
