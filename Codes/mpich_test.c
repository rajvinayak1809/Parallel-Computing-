#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 20

int main(int argc, char *argv[]) {
    MPI_File fh;  // FILE
    MPI_Offset fo;
    int myrank, nprocs;
    int buf[BUFSIZE], rbuf[BUFSIZE];
    char filename[] = "testfileIO";
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Initialize buffer
    for (int i = 0; i < BUFSIZE; i++)
        buf[i] = i;
    
    // File open, fh: individual file pointer
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

    // Compute file offset for each process
    fo = (MPI_Offset) myrank * BUFSIZE * sizeof(int);

    // File write using explicit offset (independent I/O)
    MPI_File_write_at(fh, fo, buf, BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);

    // Ensure data is written before reading
    MPI_Barrier(MPI_COMM_WORLD);

    // File read using explicit offset (independent I/O)
    MPI_File_read_at(fh, fo, rbuf, BUFSIZE, MPI_INT, &status);

    // Close file
    MPI_File_close(&fh);

    if(myrank){
        // Verify correctness
        for (int i = 0; i < BUFSIZE; i++) {
            printf("%d : rbuf[%d] = %d\n",i,i,rbuf[i]);
            if (buf[i] != rbuf[i])
                printf("Process %d: Mismatch [%d] %d %d\n", myrank, i, buf[i], rbuf[i]);
        }
    }
    
    MPI_Finalize();
    return 0;
}