#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

typedef struct {
    int x, y, z;
} Point;

Point pos_to_coords(int i, int X, int Y, int Z) {
    Point coords;
    coords.z = i / (X * Y);
    int remainder = i % (X * Y);
    coords.y = remainder / X;
    coords.x = remainder % X;
    return coords;
}

int coords_to_pos(Point coords, int X, int Y, int Z) {
    return coords.z * X * Y + coords.y * X + coords.x;
}

int read_input_file(const char* filename, double**** data, int NC, int NX, int NY, int NZ, int PX, int PY, int PZ) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file");
        return -1;
    }

    printf("Process 0: Started reading input file: %s\n", filename);

    int pos_per_rank = (NX * NY * NZ) / (PX * PY * PZ);
    int size = PX * PY * PZ;
    int NX_per_process = NX / PX;
    int NY_per_process = NY / PY;
    int NZ_per_process = NZ / PZ;

    *data = (double***)malloc(size * sizeof(double**));
    for(int i = 0; i < size; i++) {
        (*data)[i] = (double**)malloc((pos_per_rank) * sizeof(double*));
    }

    int* idx = (int*)malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        idx[i] = 0;
    }

    for (int i = 0; i < NX * NY * NZ; i++) {
        double* values = (double*)malloc(NC * sizeof(double));
        for (int c = 0; c < NC; c++) {
            if (fscanf(file, "%lf", &values[c]) != 1) {
                printf("Error reading file");
                return -1;
            }
        }

        Point data_coords = pos_to_coords(i, NX, NY, NZ);
        Point rank_coords = {
            data_coords.x / NX_per_process,
            data_coords.y / NY_per_process,
            data_coords.z / NZ_per_process
        };
        int rank = coords_to_pos(rank_coords, PX, PY, PZ);
        (*data)[rank][idx[rank]++] = values;
    }

    fclose(file);
    printf("Process 0: Completed reading file.\n");

    return 0;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 10) {
        if (rank == 0) {
            printf("Usage: %s <input_filename> <PX> <PY> <PZ> <NX> <NY> <NZ> <NC> <output_filename>\n", argv[0]);
        }
        MPI_Finalize();
        return -1;
    }

    char *input_filename = argv[1];
    int PX = atoi(argv[2]);
    int PY = atoi(argv[3]);
    int PZ = atoi(argv[4]);
    int NX = atoi(argv[5]);
    int NY = atoi(argv[6]);
    int NZ = atoi(argv[7]);
    int NC = atoi(argv[8]);
    char *output_filename = argv[9];
    

    if (PX * PY * PZ != size) {
        if (rank == 0) {
            fprintf(stderr, "Error: PX*PY*PZ must equal the number of processes\n");
        }
        MPI_Finalize();
        return -1;
    }

    double time1, time2, time3, time4;

    time1 = MPI_Wtime();
    
    int pos_per_rank = (NX *
    double** local_data = (double**)malloc(pos_per_rank * sizeof(double*));
    for(int i = 0; i < pos_per_rank; i++) {
        local_data[i] = (double*)malloc(NC * sizeof(double));
    }
    
    if (rank == 0) {
        double*** global_data = NULL;
        if (read_input_file(input_filename, &global_data, NC, NX, NY, NZ, PX, PY, PZ) != 0) {
            MPI_Finalize();
            return -1;
        }

        for (int i = 0; i < size; i++) {
            if (i == 0) {
                for (int j = 0; j < pos_per_rank; j++) {
                    memcpy(local_data[j], global_data[0][j], NC * sizeof(double));
                }
            } else {
                for(int j = 0; j < pos_per_rank; j++) {
                    MPI_Send(global_data[i][j], NC, MPI_DOUBLE, i, i*pos_per_rank + j, MPI_COMM_WORLD);
                }
            }
        }

        for (int i = 0; i < size; i++) {
            for(int j = 0; j < pos_per_rank; j++) {
                free(global_data[i][j]);
            }
            free(global_data[i]);
        }
        free(global_data);
    } else {
        for(int i = 0; i < pos_per_rank; i++) {
            MPI_Recv(local_data[i], NC, MPI_DOUBLE, 0, rank*pos_per_rank + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        printf("Process %d: Data received. Processing...\n", rank);
    }

    time2 = MPI_Wtime();


    // TODO: complete the remaning parts of the code


    free(local_data);

    MPI_Finalize();
    return 0;
}
