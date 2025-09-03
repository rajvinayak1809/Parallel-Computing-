#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <float.h>

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
    free(idx);

    return 0;
}

// Function to check if a point is a local minimum or maximum
int is_local_min_max(double* values, int NC, int local_pos, int NX, int NY, int NZ, 
                     int x, int y, int z, double* neighbors_values) {
    int local_min = 1, local_max = 1;
    
    for (int t = 0; t < NC; t++) {
        double current_val = values[t];
        
        // Check 26 neighboring points
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    
                    int nx = x + dx, ny = y + dy, nz = z + dz;
                    
                    // Skip if neighbor is out of bounds
                    if (nx < 0 || nx >= NX || ny < 0 || ny >= NY || nz < 0 || nz >= NZ) continue;
                    
                    int neighbor_pos = nz * NX * NY + ny * NX + nx;
                    
                    // Use neighbors_values for comparison
                    double neighbor_val = neighbors_values[t];
                    
                    if (current_val >= neighbor_val) local_min = 0;
                    if (current_val <= neighbor_val) local_max = 0;
                }
            }
        }
    }
    
    // Special case for points at the domain boundaries
    if (x == 0 || x == NX-1 || y == 0 || y == NY-1 || z == 0 || z == NZ-1) {
        local_min = local_max = 0;
    }
    
    return local_min ? 1 : (local_max ? 2 : 0);
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
    
    int pos_per_rank = (NX * NY * NZ) / (PX * PY * PZ);
    int NX_per_process = NX / PX;
    int NY_per_process = NY / PY;
    int NZ_per_process = NZ / PZ;

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
                local_data = global_data[0];
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
    }

    time2 = MPI_Wtime();

    // Process local domain
    int* local_min_max_counts = (int*)calloc(NC * 2, sizeof(int));
    double* local_global_extremes = (double*)malloc(NC * 2 * sizeof(double));
    
    for (int t = 0; t < NC; t++) {
        local_global_extremes[t*2] = DBL_MAX;     // Local minimum
        local_global_extremes[t*2 + 1] = -DBL_MAX; // Local maximum
    }

    // Distribute ghost cell values from neighboring processes
    // Note: This is a simplified approach and would need more complex communication in a full solution
    double* neighbors_values = (double*)malloc(NC * sizeof(double));

    Point local_rank_coords = pos_to_coords(rank, PX, PY, PZ);
    
    for (int local_pos = 0; local_pos < pos_per_rank; local_pos++) {
        // Calculate global coordinates for this point
        int local_z = local_pos / (NX_per_process * NY_per_process);
        int remainder = local_pos % (NX_per_process * NY_per_process);
        int local_y = remainder / NX_per_process;
        int local_x = remainder % NX_per_process;
        
        int global_x = local_rank_coords.x * NX_per_process + local_x;
        int global_y = local_rank_coords.y * NY_per_process + local_y;
        int global_z = local_rank_coords.z * NZ_per_process + local_z;

        // Temporary placeholder for neighbors values (would require inter-process communication)
        memcpy(neighbors_values, local_data[local_pos], NC * sizeof(double));

        int min_max_status = is_local_min_max(local_data[local_pos], NC, local_pos, NX, NY, NZ, 
                                               global_x, global_y, global_z, neighbors_values);

        for (int t = 0; t < NC; t++) {
            // Update local min/max
            if (local_data[local_pos][t] < local_global_extremes[t*2]) 
                local_global_extremes[t*2] = local_data[local_pos][t];
            if (local_data[local_pos][t] > local_global_extremes[t*2 + 1]) 
                local_global_extremes[t*2 + 1] = local_data[local_pos][t];

            // Count local minima and maxima
            if (min_max_status == 1) local_min_max_counts[t*2]++;
            else if (min_max_status == 2) local_min_max_counts[t*2 + 1]++;
        }
    }

    // Collect results from all processes
    int* global_min_max_counts = NULL;
    double* global_extremes = NULL;

    if (rank == 0) {
        global_min_max_counts = (int*)calloc(NC * 2, sizeof(int));
        global_extremes = (double*)malloc(NC * 2 * sizeof(double));
        
        for (int t = 0; t < NC; t++) {
            global_extremes[t*2] = DBL_MAX;
            global_extremes[t*2 + 1] = -DBL_MAX;
        }
    }

    // Reduce min/max counts
    MPI_Reduce(local_min_max_counts, global_min_max_counts, NC * 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
    // Reduce global extremes
    MPI_Reduce(local_global_extremes, global_extremes, NC * 2, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_global_extremes, global_extremes + 1, NC * 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    time3 = MPI_Wtime();

    // Write output from rank 0
    if (rank == 0) {
        // Timing information
        time4 = MPI_Wtime();
        double read_time = time2 - time1;
        double comp_time = time3 - time2;
        double total_time = time4 - time1;
        
        double max_time1, max_time2, max_time3;
        
        MPI_Reduce(&read_time, &max_time1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&comp_time, &max_time2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_time, &max_time3, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        FILE* output_file = fopen(output_filename, "w");
        
        // Line 1: Local minima and maxima
        for (int t = 0; t < NC; t++) {
            fprintf(output_file, "(%d,%d)%s", 
                    global_min_max_counts[t*2], 
                    global_min_max_counts[t*2 + 1],
                    t < NC-1 ? ", " : "\n");
        }
        
        // Line 2: Global extremes
        for (int t = 0; t < NC; t++) {
            fprintf(output_file, "(%f,%f)%s", 
                    global_extremes[t*2], 
                    global_extremes[t*2 + 1],
                    t < NC-1 ? ", " : "\n");
        }
        
        // Line 3: Timing information
        fprintf(output_file, "%f, %f, %f\n", max_time1, max_time2, max_time3);
        
        fclose(output_file);
    }

    // Cleanup
    free(local_data);
    free(local_min_max_counts);
    free(local_global_extremes);
    free(neighbors_values);

    if (rank == 0) {
        free(global_min_max_counts);
        free(global_extremes);
    }

    MPI_Finalize();
    return 0;
}