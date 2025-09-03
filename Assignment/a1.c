#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <float.h>

void find_local_extrema(double *data, int local_size, int time_steps, int *local_min_count, int *local_max_count, double *local_min, double *local_max) {
    *local_min = DBL_MAX;
    *local_max = -DBL_MAX;
    *local_min_count = 0;
    *local_max_count = 0;
    
    for (int t = 0; t < time_steps; t++) {
        for (int i = 1; i < local_size - 1; i++) { 
            double prev = data[(i - 1) * time_steps + t];
            double curr = data[i * time_steps + t];
            double next = data[(i + 1) * time_steps + t];
            
            if (curr < prev && curr < next) (*local_min_count)++;
            if (curr > prev && curr > next) (*local_max_count)++;
            if (curr < *local_min) *local_min = curr;
            if (curr > *local_max) *local_max = curr;
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 10) {
        if (rank == 0) printf("Usage: %s <input_file> <PX> <PY> <PZ> <NX> <NY> <NZ> <NC> <output_file>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }
    
    char *input_file = argv[1];
    int PX = atoi(argv[2]), PY = atoi(argv[3]), PZ = atoi(argv[4]);
    int NX = atoi(argv[5]), NY = atoi(argv[6]), NZ = atoi(argv[7]), NC = atoi(argv[8]);
    char *output_file = argv[9];
    
    int total_points = NX * NY * NZ;
    int local_size = total_points / size;
    double *data = NULL;
    double *local_data = (double *)malloc(local_size * NC * sizeof(double));
    
    if (rank == 0) {
        data = (double *)malloc(total_points * NC * sizeof(double));
        FILE *file = fopen(input_file, "r");
        if (!file) {
            printf("Error opening file\n");
            MPI_Finalize();
            return -1;
        }
        
        for (int i = 0; i < total_points * NC; i++) {
            fscanf(file, "%lf", &data[i]);
        }
        fclose(file);
        printf("Process 0: Completed reading %d rows from input file.\n", total_points);
    }
    
    MPI_Scatter(data, local_size * NC, MPI_DOUBLE, local_data, local_size * NC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) free(data);
    printf("Process %d: Received %d data points.\n", rank, local_size);
    
    int local_min_count, local_max_count;
    double local_min, local_max;
    find_local_extrema(local_data, local_size, NC, &local_min_count, &local_max_count, &local_min, &local_max);
    
    int global_min_count, global_max_count;
    double global_min, global_max;
    
    MPI_Reduce(&local_min_count, &global_min_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_max_count, &global_max_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        FILE *out = fopen(output_file, "w");
        fprintf(out, "Global Min: %lf\n", global_min);
        fprintf(out, "Global Max: %lf\n", global_max);
        fprintf(out, "Total Local Minima Count: %d\n", global_min_count);
        fprintf(out, "Total Local Maxima Count: %d\n", global_max_count);
        fclose(out);
    }
    
    free(local_data);
    MPI_Finalize();
    return 0;
}
