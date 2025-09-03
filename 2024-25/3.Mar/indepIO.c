// Set explicit file offset

#include "mpi.h" 
#include <stdio.h> 
#include <string.h> 
#include <stdlib.h> 

int main(int argc, char *argv[]) {    

  int i, myrank, nprocs, numbytes;
  int BUFSIZE = atoi(argv[1]);
  int buf[BUFSIZE], rbuf[BUFSIZE], tsize;
  char filename[128];
  MPI_Status status;
  double tDiff;
  double MB=1048576.0;

  MPI_File fh;  // FILE

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Assuming every process writes 'numbytes' bytes (BUFSIZE integers) 
  numbytes = BUFSIZE*sizeof(int);
  double fileSize=nprocs*numbytes/MB;

  for (int i=0; i<BUFSIZE ; i++)
    buf[i]=i;

  strcpy(filename, "testfileIO");

  // File open, fh: individual file pointer
  MPI_File_open (MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

  MPI_Offset fo = (MPI_Offset) myrank*BUFSIZE*sizeof(int);

  double tStart = MPI_Wtime();
  // File write using explicit offset (independent I/O)
  MPI_File_write_at (fh, fo, buf, BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);	
  double tEnd = MPI_Wtime() - tStart;
  MPI_Reduce(&tEnd, &tDiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) printf ("Time to write %lf MB file = %lf seconds | Write bandwidth = %lf\n", tDiff, fileSize, fileSize/tDiff);

  tStart = MPI_Wtime();
  // File read using explicit offset (independent I/O)
  MPI_File_read_at (fh, fo, rbuf, BUFSIZE, MPI_INT, &status);	
  tEnd = MPI_Wtime() - tStart;
  MPI_Reduce(&tEnd, &tDiff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) printf ("Time to read %lf MB file = %lf seconds | Read bandwidth = %lf\n", tDiff, fileSize, fileSize/tDiff);

  MPI_File_close (&fh);		//fclose

  for (i=0; i<BUFSIZE ; i++)
    if (buf[i] != rbuf[i]) printf ("Mismatch [%d] %d %d\n", i, buf[i], rbuf[i]);

  MPI_Finalize();    
  return 0; 

} 

