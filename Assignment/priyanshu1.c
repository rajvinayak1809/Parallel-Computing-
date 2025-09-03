#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <math.h>

typedef struct
{
    int x, y, z;
} Point;

Point pos_to_coords(int i, int X, int Y, int Z)
{
    Point coords;
    coords.z = i / (X * Y);
    int remainder = i % (X * Y);
    coords.y = remainder / X;
    coords.x = remainder % X;
    return coords;
}

int coords_to_pos(Point coords, int X, int Y, int Z)
{
    return coords.z * X * Y + coords.y * X + coords.x;
}

int read_input_file(const char *filename, double ***data, int NC, int NX, int NY, int NZ, int PX, int PY, int PZ)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error opening file");
        return -1;
    }

    int pos_per_rank = (NX * NY * NZ) / (PX * PY * PZ);
    int size = PX * PY * PZ;
    int NX_per_process = NX / PX;
    int NY_per_process = NY / PY;
    int NZ_per_process = NZ / PZ;

    *data = (double **)malloc(NC * sizeof(double *));
    for (int i = 0; i < NC; i++)
    {
        (*data)[i] = (double *)malloc(NX * NY * NZ * sizeof(double));
    }

    int *idx = (int *)malloc(size * sizeof(int));
    for (int i = 0; i < size; i++)
    {
        idx[i] = 0;
    }

    for (int i = 0; i < NX * NY * NZ; i++)
    {
        Point data_coords = pos_to_coords(i, NX, NY, NZ);
        Point rank_coords = {
            data_coords.x / NX_per_process,
            data_coords.y / NY_per_process,
            data_coords.z / NZ_per_process};
        int rank = coords_to_pos(rank_coords, PX, PY, PZ);

        for (int c = 0; c < NC; c++)
        {
            double val;
            if (fscanf(file, "%lf", &val) != 1)
            {
                printf("Error reading file");
                return -1;
            }
            (*data)[c][rank * pos_per_rank + idx[rank]] = val;
        }

        idx[rank]++;
    }

    fclose(file);

    return 0;
}

bool isLocalMinimum(double *local_data, int local_index, int NC, int NX, int NY, int NZ, int PX, int PY, int PZ, double *X_plus_ones_layer_recv, double *X_minus_plus_ones_layer_recv, Point rankCoordinates)
{
    Point local_coordinates = pos_to_coords(local_index, NX / PX, NY / PY, NZ / PZ);

    if (local_coordinates.x == 0 && rankCoordinates.x > 0)
    {
        int neighbour_index = (NY / PY) * local_coordinates.z + local_coordinates.y;
        if (X_minus_plus_ones_layer_recv[neighbour_index] < local_data[local_index])
            return false;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 10)
    {
        if (rank == 0)
        {
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

    if (PX * PY * PZ != size)
    {
        if (rank == 0)
        {
            // printing size and px*py*pz
            fprintf(stderr, "Error: Number of processes (%d) does not match the product of PX, PY, and PZ (%d)\n", size, PX * PY * PZ);
            fprintf(stderr, "Error: PX*PY*PZ must equal the number of processes\n");
        }
        MPI_Finalize();
        return -1;
    }

    double time1, time2, time3, time4;

    time1 = MPI_Wtime();

    int elements_per_rank = (NX * NY * NZ) / (PX * PY * PZ);
    int NX_per_process = NX / PX;
    int NY_per_process = NY / PY;
    int NZ_per_process = NZ / PZ;

    double **local_data = (double **)malloc(NC * sizeof(double *));
    for (int i = 0; i < NC; i++)
    {
        local_data[i] = (double *)malloc(elements_per_rank * sizeof(double));
    }

    if (rank == 0)
    {
        double **global_data = NULL;
        if (read_input_file(input_filename, &global_data, NC, NX, NY, NZ, PX, PY, PZ) != 0)
        {
            MPI_Finalize();
            return -1;
        }

        for (int i = 0; i < NC; i++)
        {
            MPI_Scatter(global_data[i], elements_per_rank, MPI_DOUBLE, local_data[i], elements_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            free(global_data[i]);
        }

        free(global_data);
    }
    else
    {
        for (int i = 0; i < NC; i++)
        {
            MPI_Scatter(NULL, elements_per_rank, MPI_DOUBLE, local_data[i], elements_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // printing the local data received by each process
    printf("Rank %d received data:\n", rank);
    for (int i = 0; i < NC; i++)
    {
        printf("Data for testcase %d:\n", i);
        for (int j = 0; j < elements_per_rank; j++)
        {
            printf("%lf ", local_data[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

    MPI_Barrier(MPI_COMM_WORLD);

    time2 = MPI_Wtime();

    double **X_plus_ones_layer_send = (double **)malloc(NC * sizeof(double *));
    double **X_minus_ones_layer_send = (double **)malloc(NC * sizeof(double *));
    double **Y_plus_ones_layer_send = (double **)malloc(NC * sizeof(double *));
    double **Y_minus_ones_layer_send = (double **)malloc(NC * sizeof(double *));
    double **Z_plus_ones_layer_send = (double **)malloc(NC * sizeof(double *));
    double **Z_minus_ones_layer_send = (double **)malloc(NC * sizeof(double *));

    double **X_plus_ones_layer_recv = (double **)malloc(NC * sizeof(double *));
    double **X_minus_ones_layer_recv = (double **)malloc(NC * sizeof(double *));
    double **Y_plus_ones_layer_recv = (double **)malloc(NC * sizeof(double *));
    double **Y_minus_ones_layer_recv = (double **)malloc(NC * sizeof(double *));
    double **Z_plus_ones_layer_recv = (double **)malloc(NC * sizeof(double *));
    double **Z_minus_ones_layer_recv = (double **)malloc(NC * sizeof(double *));

    // Allocate memory for the layers
    for (int i = 0; i < NC; i++)
    {
        X_plus_ones_layer_send[i] = (double *)malloc(elements_per_rank / NX_per_process * sizeof(double));
        X_minus_ones_layer_send[i] = (double *)malloc(elements_per_rank / NX_per_process * sizeof(double));
        Y_plus_ones_layer_send[i] = (double *)malloc(elements_per_rank / NY_per_process * sizeof(double));
        Y_minus_ones_layer_send[i] = (double *)malloc(elements_per_rank / NY_per_process * sizeof(double));
        Z_plus_ones_layer_send[i] = (double *)malloc(elements_per_rank / NZ_per_process * sizeof(double));
        Z_minus_ones_layer_send[i] = (double *)malloc(elements_per_rank / NZ_per_process * sizeof(double));

        X_plus_ones_layer_recv[i] = (double *)malloc(elements_per_rank / NX_per_process * sizeof(double));
        X_minus_ones_layer_recv[i] = (double *)malloc(elements_per_rank / NX_per_process * sizeof(double));
        Y_plus_ones_layer_recv[i] = (double *)malloc(elements_per_rank / NY_per_process * sizeof(double));
        Y_minus_ones_layer_recv[i] = (double *)malloc(elements_per_rank / NY_per_process * sizeof(double));
        Z_plus_ones_layer_recv[i] = (double *)malloc(elements_per_rank / NZ_per_process * sizeof(double));
        Z_minus_ones_layer_recv[i] = (double *)malloc(elements_per_rank / NZ_per_process * sizeof(double));
    }

    Point rank_coordinates = pos_to_coords(rank, PX, PY, PZ);

    // communication with X-1 rank
    if (rank_coordinates.x > 0)
    {
        Point X_minus_one_coords = {rank_coordinates.x - 1, rank_coordinates.y, rank_coordinates.z};
        int X_minus_one_rank = coords_to_pos(X_minus_one_coords, PX, PY, PZ);
        for (int t = 0; t < NC; t++)
        {
            int index = 0;
            for (int i = 0; i < elements_per_rank; i += NX_per_process)
            {
                X_minus_ones_layer_send[t][index++] = local_data[t][i];
            }
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Send(X_minus_ones_layer_send[t], elements_per_rank / NX_per_process, MPI_DOUBLE, X_minus_one_rank, 0, MPI_COMM_WORLD);
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Recv(X_minus_ones_layer_recv[t], elements_per_rank / NX_per_process, MPI_DOUBLE, X_minus_one_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // printing the received layer
        for (int i = 0; i < NC; i++)
        {
            printf("Rank %d received X-1 layer from rank %d for testcase %d:\n", rank, X_minus_one_rank, i);
            for (int j = 0; j < elements_per_rank / NX_per_process; j++)
            {
                printf("%lf ", X_minus_ones_layer_recv[i][j]);
            }
            printf("\n");
        }
    }

    // communication with X+1 rank
    if (rank_coordinates.x < PX - 1)
    {

        Point X_plus_one_coords = {rank_coordinates.x + 1, rank_coordinates.y, rank_coordinates.z};
        int X_plus_one_rank = coords_to_pos(X_plus_one_coords, PX, PY, PZ);

        for (int t = 0; t < NC; t++)
        {
            MPI_Recv(X_plus_ones_layer_recv[t], elements_per_rank / NX_per_process, MPI_DOUBLE, X_plus_one_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int t = 0; t < NC; t++)
        {
            int index = 0;
            for (int i = NX_per_process - 1; i < elements_per_rank; i += NX_per_process)
            {
                X_plus_ones_layer_send[t][index++] = local_data[t][i];
            }
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Send(X_plus_ones_layer_send[t], elements_per_rank / NX_per_process, MPI_DOUBLE, X_plus_one_rank, 0, MPI_COMM_WORLD);
        }

        // printing the received layer
        for (int i = 0; i < NC; i++)
        {
            printf("Rank %d received X+1 layer from rank %d for testcase %d:\n", rank, X_plus_one_rank, i);
            for (int j = 0; j < elements_per_rank / NX_per_process; j++)
            {
                printf("%lf ", X_plus_ones_layer_recv[i][j]);
            }
            printf("\n");
        }
    }

    // communication with Y-1 rank
    if (rank_coordinates.y > 0)
    {
        Point Y_minus_one_coords = {rank_coordinates.x, rank_coordinates.y - 1, rank_coordinates.z};
        int Y_minus_one_rank = coords_to_pos(Y_minus_one_coords, PX, PY, PZ);

        for (int t = 0; t < NC; t++)
        {
            int index = 0;
            for (int z = 0; z < NZ_per_process; z++)
            {
                for (int x = 0; x < NX_per_process; x++)
                {
                    Y_minus_ones_layer_send[t][index++] = local_data[t][z * NX_per_process * NY_per_process + x];
                }
            }
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Send(Y_minus_ones_layer_send[t], elements_per_rank / NY_per_process, MPI_DOUBLE, Y_minus_one_rank, 0, MPI_COMM_WORLD);
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Recv(Y_minus_ones_layer_recv[t], elements_per_rank / NY_per_process, MPI_DOUBLE, Y_minus_one_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // printing the received layer
        for (int i = 0; i < NC; i++)
        {
            printf("Rank %d received Y-1 layer from rank %d for testcase %d:\n", rank, Y_minus_one_rank, i);
            for (int j = 0; j < elements_per_rank / NY_per_process; j++)
            {
                printf("%lf ", Y_minus_ones_layer_recv[i][j]);
            }
            printf("\n");
        }
    }

    // communication with Y+1 rank
    if (rank_coordinates.y < PY - 1)
    {
        Point Y_plus_one_coords = {rank_coordinates.x, rank_coordinates.y + 1, rank_coordinates.z};
        int Y_plus_one_rank = coords_to_pos(Y_plus_one_coords, PX, PY, PZ);

        for (int t = 0; t < NC; t++)
        {
            MPI_Recv(Y_plus_ones_layer_recv[t], elements_per_rank / NY_per_process, MPI_DOUBLE, Y_plus_one_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int t = 0; t < NC; t++)
        {
            int index = 0;
            for (int z = 0; z < NZ_per_process; z++)
            {
                for (int x = 0; x < NX_per_process; x++)
                {
                    Y_plus_ones_layer_send[t][index++] = local_data[t][z * NX_per_process * NY_per_process + (NY_per_process - 1) * NX_per_process + x];
                }
            }
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Send(Y_plus_ones_layer_send[t], elements_per_rank / NY_per_process, MPI_DOUBLE, Y_plus_one_rank, 0, MPI_COMM_WORLD);
        }
        // printing the received layer
        for (int i = 0; i < NC; i++)
        {
            printf("Rank %d received Y+1 layer from rank %d for testcase %d:\n", rank, Y_plus_one_rank, i);
            for (int j = 0; j < elements_per_rank / NY_per_process; j++)
            {
                printf("%lf ", Y_plus_ones_layer_recv[i][j]);
            }
            printf("\n");
        }
    }

    // communication with Z-1 rank
    if (rank_coordinates.z > 0)
    {
        Point Z_minus_one_coords = {rank_coordinates.x, rank_coordinates.y, rank_coordinates.z - 1};
        int Z_minus_one_rank = coords_to_pos(Z_minus_one_coords, PX, PY, PZ);

        for (int t = 0; t < NC; t++)
        {
            int index = 0;
            for (int i = 0; i < NX_per_process * NY_per_process; i++)
            {
                Z_minus_ones_layer_send[t][index++] = local_data[t][i];
            }
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Send(Z_minus_ones_layer_send[t], elements_per_rank / NZ_per_process, MPI_DOUBLE, Z_minus_one_rank, 0, MPI_COMM_WORLD);
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Recv(Z_minus_ones_layer_recv[t], elements_per_rank / NZ_per_process, MPI_DOUBLE, Z_minus_one_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // printing the received layer
        for (int i = 0; i < NC; i++)
        {
            printf("Rank %d received Z-1 layer from rank %d for testcase %d:\n", rank, Z_minus_one_rank, i);
            for (int j = 0; j < elements_per_rank / NZ_per_process; j++)
            {
                printf("%lf ", Z_minus_ones_layer_recv[i][j]);
            }
            printf("\n");
        }
    }

    // communication with Z+1 rank
    if (rank_coordinates.z < PZ - 1)
    {
        Point Z_plus_one_coords = {rank_coordinates.x, rank_coordinates.y, rank_coordinates.z + 1};
        int Z_plus_one_rank = coords_to_pos(Z_plus_one_coords, PX, PY, PZ);

        for (int t = 0; t < NC; t++)
        {
            MPI_Recv(Z_plus_ones_layer_recv[t], elements_per_rank / NZ_per_process, MPI_DOUBLE, Z_plus_one_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int t = 0; t < NC; t++)
        {
            int index = 0;
            for (int i = 0; i < NX_per_process * NY_per_process; i++)
            {
                Z_plus_ones_layer_send[t][index++] = local_data[t][i + (NZ_per_process - 1) * NX_per_process * NY_per_process];
            }
        }

        for (int t = 0; t < NC; t++)
        {
            MPI_Send(Z_plus_ones_layer_send[t], elements_per_rank / NZ_per_process, MPI_DOUBLE, Z_plus_one_rank, 0, MPI_COMM_WORLD);
        }

        // printing the received layer
        for (int i = 0; i < NC; i++)
        {
            printf("Rank %d received Z+1 layer from rank %d for testcase %d:\n", rank, Z_plus_one_rank, i);
            for (int j = 0; j < elements_per_rank / NZ_per_process; j++)
            {
                printf("%lf ", Z_plus_ones_layer_recv[i][j]);
            }
            printf("\n");
        }
    }
    time3 = MPI_Wtime();

    // Checking for the local minimum and maximums
    int local_minimaCount[NC], local_maximaCount[NC];
    double local_maximum[NC], local_minimum[NC];
    for (int i = 0; i < NC; i++)
    {
        local_maximum[i] = local_data[i][0];
        local_minimum[i] = local_data[i][0];
        local_minimaCount[i] = 0;
        local_maximaCount[i] = 0;
    }

    for (int t = 0; t < NC; t++)
    {
        for (int i = 0; i < elements_per_rank; i++)
        {

            // if(isLocalMinimum()){
            //     local_minimaCount[t]++;
            // }
            // if(isLocalMaximum()){
            //     local_maximaCount[t]++;
            // }

            local_maximum[t] = fmax(local_maximum[t], local_data[t][i]);
            local_minimum[t] = fmin(local_minimum[t], local_data[t][i]);
        }
    }

    /*
    Output: The output file should contain three lines (no extra lines to be printed). The number of 2-tuples in lines 1 and 2
    depends on the number of time steps.
    Line 1: (count of local minima, count of local maxima), ...
    Line 2: (global minimum, global maximum), ...
    Line 3: Read Time, Main code Time, Total time // a pair of integer values per time step
    // a pair of float values per time step
    // three doubles (across all time steps)
     */

    // writing to the output file
    if (!rank)
    {
        FILE *output_file = fopen(output_filename, "w");
        if (output_file == NULL)
        {
            printf("Error opening output file");
            return -1;
        }
        for (int i = 0; i < NC; i++)
        {
            fprintf(output_file, "(%d, %d), ", local_minimaCount[i], local_maximaCount[i]);
        }
        fprintf(output_file, "\n");
        for (int i = 0; i < NC; i++)
        {
            fprintf(output_file, "(%lf, %lf), ", local_minimum[i], local_maximum[i]);
        }
        fprintf(output_file, "\n");
        fprintf(output_file, "%lf %lf %lf\n", time2 - time1, time3 - time2, time3 - time1);
        fclose(output_file);
        printf("Output written to %s\n", output_filename);
    }

    free(local_data);

    MPI_Finalize();
    return 0;
}