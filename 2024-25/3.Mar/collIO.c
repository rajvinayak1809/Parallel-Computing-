#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int myrank, commsize;
int size, numbytes;
int *sdbuf;
double tStart, tEnd, tDiff;
double MB=1048576.0;
MPI_Datatype etype, filetype;

void prnerror (int error, char *string)
{
    fprintf(stderr, "Error %d in %s\n", error, string);
    MPI_Finalize();
    exit(-1);
}

void file_write () {

  MPI_File fileHandle;
  MPI_Status status;
  char fileNameFS[16]; // = "dummy";
  sprintf(fileNameFS, "dummy-%d", commsize);
  int mode = MPI_MODE_CREATE | MPI_MODE_RDWR;

  MPI_File_open (MPI_COMM_WORLD, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
  MPI_Offset disp = (MPI_Offset) myrank*numbytes;

  etype = MPI_INT;
  filetype = MPI_INT;
  MPI_File_set_view(fileHandle, disp, etype, filetype, "native", MPI_INFO_NULL);
  tStart = MPI_Wtime();
  int result = MPI_File_write_all (fileHandle, sdbuf, size, MPI_INT, &status);
  tEnd = MPI_Wtime() - tStart;
  MPI_Reduce(&tEnd, &tDiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_File_close (&fileHandle);
  if (result != MPI_SUCCESS) 
       prnerror (result, "MPI_File_write Error:");

  double filesize=commsize*numbytes/MB;
  if (myrank == 0)
    printf("Time to write %8.3lf MB = %lf seconds WB = %9.3lf MB/s\n", filesize, tDiff, filesize/tDiff);

}

void file_read () {

  MPI_File fileHandle;
  MPI_Status status;
  char fileNameFS[16]; // = "dummy";
  sprintf(fileNameFS, "dummy-%d", commsize);
  int mode = MPI_MODE_RDWR ; 

  MPI_Datatype etype, filetype, twoints;
  MPI_File_open (MPI_COMM_WORLD, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
  MPI_Offset disp = (MPI_Offset)myrank*numbytes;

  etype = MPI_INT;
  MPI_File_set_view(fileHandle, disp, etype, filetype, "native", MPI_INFO_NULL);
  tStart = MPI_Wtime();
  int result = MPI_File_read_all (fileHandle, sdbuf, size, MPI_INT, &status);
  tEnd = MPI_Wtime() - tStart;
  MPI_Reduce(&tEnd, &tDiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_File_close (&fileHandle);
  if (result != MPI_SUCCESS) 
       prnerror (result, "MPI_Read_write Error:");

  double filesize=commsize*numbytes/MB;
  if (myrank == 0)
    printf("Time to read %8.3lf MB = %lf seconds RB = %9.3lf MB/s\n", filesize, tDiff, filesize/tDiff);

}

int main(int argc, char *argv[]) {

  long int i;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  MPI_Comm_size (MPI_COMM_WORLD, &commsize);

  size = atoi(argv[1]);
  numbytes = size * sizeof(int);

  posix_memalign((void **)&sdbuf, 64, numbytes);
  if (!sdbuf) {
    printf("%d: Error allocating memory %d bytes failed\n", myrank, numbytes);
    fflush(stdout);
    exit(1);
  }

  //Random initialization
  for (i=0; i<size; i++) 
    sdbuf[i] = i+5 * 3 * myrank;

  file_write();
  file_read();

  MPI_Finalize();

}

