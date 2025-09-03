#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <float.h>
#include <math.h>

// Function to read input file
void readInputFile(const char *filename, float *data, int nx, int ny, int nz, int nc)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Error opening file: %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < nx * ny * nz; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            if (fscanf(file, "%f", &data[i * nc + j]) != 1)
            {
                fprintf(stderr, "Error reading data at point %d, timestep %d\n", i, j);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    fclose(file);
}

// Function to determine if a point is a local minimum
int isLocalMinimum(float *data, int x, int y, int z, int t, int nx, int ny, int nz, int nc,
                   int local_nx, int local_ny, int local_nz,
                   int offset_x, int offset_y, int offset_z,
                   float *ghost_layers)
{

    float center_value = data[(z * local_ny * local_nx + y * local_nx + x) * nc + t];
    int is_minimum = 1;

    // Check 6 neighbors (3D von Neumann neighborhood)
    // Directions: -x, +x, -y, +y, -z, +z
    int dx[6] = {-1, 1, 0, 0, 0, 0};
    int dy[6] = {0, 0, -1, 1, 0, 0};
    int dz[6] = {0, 0, 0, 0, -1, 1};

    for (int i = 0; i < 6; i++)
    {
        int nx_check = x + dx[i];
        int ny_check = y + dy[i];
        int nz_check = z + dz[i];

        // Skip check if neighbor is outside the global domain
        if ((offset_x + nx_check < 0) || (offset_x + nx_check >= nx) ||
            (offset_y + ny_check < 0) || (offset_y + ny_check >= ny) ||
            (offset_z + nz_check < 0) || (offset_z + nz_check >= nz))
        {
            continue;
        }

        float neighbor_value;

        // If the neighbor is within the local domain
        if (nx_check >= 0 && nx_check < local_nx &&
            ny_check >= 0 && ny_check < local_ny &&
            nz_check >= 0 && nz_check < local_nz)
        {
            neighbor_value = data[(nz_check * local_ny * local_nx + ny_check * local_nx + nx_check) * nc + t];
        }
        else
        {
            // Use the ghost layer data
            // Calculate index in ghost_layers
            int ghost_index = -1;

            // Indices for ghost layers
            if (nx_check < 0)
            {
                ghost_index = 0;
            }
            else if (nx_check >= local_nx)
            {
                ghost_index = 1;
            }
            else if (ny_check < 0)
            {
                ghost_index = 2;
            }
            else if (ny_check >= local_ny)
            {
                ghost_index = 3;
            }
            else if (nz_check < 0)
            {
                ghost_index = 4;
            }
            else if (nz_check >= local_nz)
            {
                ghost_index = 5;
            }

            if (ghost_index == 0)
            { // -x face
                neighbor_value = ghost_layers[(ghost_index * local_ny * local_nz + z * local_ny + y) * nc + t];
            }
            else if (ghost_index == 1)
            { // +x face
                neighbor_value = ghost_layers[(ghost_index * local_ny * local_nz + z * local_ny + y) * nc + t];
            }
            else if (ghost_index == 2)
            { // -y face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_nz + z * local_nx + x) * nc + t];
            }
            else if (ghost_index == 3)
            { // +y face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_nz + z * local_nx + x) * nc + t];
            }
            else if (ghost_index == 4)
            { // -z face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_ny + y * local_nx + x) * nc + t];
            }
            else if (ghost_index == 5)
            { // +z face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_ny + y * local_nx + x) * nc + t];
            }
        }

        // If any neighbor has a smaller value, it's not a minimum
        if (neighbor_value < center_value)
        {
            is_minimum = 0;
            break;
        }

        // If any neighbor has the same value, it's not a strict minimum
        if (fabs(neighbor_value - center_value) < 1e-6)
        {
            is_minimum = 0;
            break;
        }
    }

    return is_minimum;
}

// Function to determine if a point is a local maximum
int isLocalMaximum(float *data, int x, int y, int z, int t, int nx, int ny, int nz, int nc,
                   int local_nx, int local_ny, int local_nz,
                   int offset_x, int offset_y, int offset_z,
                   float *ghost_layers)
{

    float center_value = data[(z * local_ny * local_nx + y * local_nx + x) * nc + t];
    int is_maximum = 1;

    // Check 6 neighbors (3D von Neumann neighborhood)
    // Directions: -x, +x, -y, +y, -z, +z
    int dx[6] = {-1, 1, 0, 0, 0, 0};
    int dy[6] = {0, 0, -1, 1, 0, 0};
    int dz[6] = {0, 0, 0, 0, -1, 1};

    for (int i = 0; i < 6; i++)
    {
        int nx_check = x + dx[i];
        int ny_check = y + dy[i];
        int nz_check = z + dz[i];

        // Skip check if neighbor is outside the global domain
        if ((offset_x + nx_check < 0) || (offset_x + nx_check >= nx) ||
            (offset_y + ny_check < 0) || (offset_y + ny_check >= ny) ||
            (offset_z + nz_check < 0) || (offset_z + nz_check >= nz))
        {
            continue;
        }

        float neighbor_value;

        // If the neighbor is within the local domain
        if (nx_check >= 0 && nx_check < local_nx &&
            ny_check >= 0 && ny_check < local_ny &&
            nz_check >= 0 && nz_check < local_nz)
        {
            neighbor_value = data[(nz_check * local_ny * local_nx + ny_check * local_nx + nx_check) * nc + t];
        }
        else
        {
            // Use the ghost layer data
            // Calculate index in ghost_layers
            int ghost_index = -1;

            // Indices for ghost layers
            if (nx_check < 0)
            {
                ghost_index = 0;
            }
            else if (nx_check >= local_nx)
            {
                ghost_index = 1;
            }
            else if (ny_check < 0)
            {
                ghost_index = 2;
            }
            else if (ny_check >= local_ny)
            {
                ghost_index = 3;
            }
            else if (nz_check < 0)
            {
                ghost_index = 4;
            }
            else if (nz_check >= local_nz)
            {
                ghost_index = 5;
            }

            if (ghost_index == 0)
            { // -x face
                neighbor_value = ghost_layers[(ghost_index * local_ny * local_nz + z * local_ny + y) * nc + t];
            }
            else if (ghost_index == 1)
            { // +x face
                neighbor_value = ghost_layers[(ghost_index * local_ny * local_nz + z * local_ny + y) * nc + t];
            }
            else if (ghost_index == 2)
            { // -y face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_nz + z * local_nx + x) * nc + t];
            }
            else if (ghost_index == 3)
            { // +y face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_nz + z * local_nx + x) * nc + t];
            }
            else if (ghost_index == 4)
            { // -z face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_ny + y * local_nx + x) * nc + t];
            }
            else if (ghost_index == 5)
            { // +z face
                neighbor_value = ghost_layers[(ghost_index * local_nx * local_ny + y * local_nx + x) * nc + t];
            }
        }

        // If any neighbor has a larger value, it's not a maximum
        if (neighbor_value > center_value)
        {
            is_maximum = 0;
            break;
        }

        // If any neighbor has the same value, it's not a strict maximum
        if (fabs(neighbor_value - center_value) < 1e-6)
        {
            is_maximum = 0;
            break;
        }
    }

    return is_maximum;
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double time1, time2, time3, time4;
    double read_time, main_time, total_time;
    double max_read_time, max_main_time, max_total_time;

    // Start timing
    time1 = MPI_Wtime();

    // Check for correct number of arguments
    if (argc != 10)
    {
        if (rank == 0)
        {
            fprintf(stderr, "Usage: %s <input_file> <PX> <PY> <PZ> <NX> <NY> <NZ> <NC> <output_file>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    // Parse command line arguments
    const char *input_file = argv[1];
    int px = atoi(argv[2]);
    int py = atoi(argv[3]);
    int pz = atoi(argv[4]);
    int nx = atoi(argv[5]);
    int ny = atoi(argv[6]);
    int nz = atoi(argv[7]);
    int nc = atoi(argv[8]);
    const char *output_file = argv[9];

    // Verify process decomposition
    if (px * py * pz != size)
    {
        if (rank == 0)
        {
            fprintf(stderr, "Error: Number of processes (%d) does not match PX*PY*PZ (%d*%d*%d=%d)\n",
                    size, px, py, pz, px * py * pz);
        }
        MPI_Finalize();
        return 1;
    }

    // Calculate process coordinates (px, py, pz)
    int coord_x = rank % px;
    int coord_y = (rank / px) % py;
    int coord_z = rank / (px * py);

    // Calculate local domain size and offsets
    int local_nx = nx / px + ((coord_x < nx % px) ? 1 : 0);
    int local_ny = ny / py + ((coord_y < ny % py) ? 1 : 0);
    int local_nz = nz / pz + ((coord_z < nz % pz) ? 1 : 0);

    // Calculate offsets
    int offset_x = (nx / px) * coord_x + ((coord_x < nx % px) ? coord_x : nx % px);
    int offset_y = (ny / py) * coord_y + ((coord_y < ny % py) ? coord_y : ny % py);
    int offset_z = (nz / pz) * coord_z + ((coord_z < nz % pz) ? coord_z : nz % pz);

    // Allocate memory for global data (only on rank 0)
    float *global_data = NULL;
    if (rank == 0)
    {
        global_data = (float *)malloc(nx * ny * nz * nc * sizeof(float));
        if (!global_data)
        {
            fprintf(stderr, "Error: Memory allocation failed for global data\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Read input file
        readInputFile(input_file, global_data, nx, ny, nz, nc);
    }

    // Allocate memory for local data
    float *local_data = (float *)malloc(local_nx * local_ny * local_nz * nc * sizeof(float));
    if (!local_data)
    {
        fprintf(stderr, "Error: Memory allocation failed for local data on rank %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Create MPI datatypes for local domain
    MPI_Datatype local_domain_type;
    int local_sizes[3] = {local_nz, local_ny, local_nx};
    int global_sizes[3] = {nz, ny, nx};
    int starts[3] = {offset_z, offset_y, offset_x};

    MPI_Type_create_subarray(3, global_sizes, local_sizes, starts, MPI_ORDER_C, MPI_FLOAT, &local_domain_type);
    MPI_Type_commit(&local_domain_type);

    // Distribute data from rank 0 to all processes
    // Need to consider all time steps for each point
    for (int t = 0; t < nc; t++)
    {
        if (rank == 0)
        {
            // Pack the data for this time step
            float *time_step_data = (float *)malloc(nx * ny * nz * sizeof(float));
            for (int i = 0; i < nx * ny * nz; i++)
            {
                time_step_data[i] = global_data[i * nc + t];
            }

            // Self copy for rank 0
            for (int z = 0; z < local_nz; z++)
            {
                for (int y = 0; y < local_ny; y++)
                {
                    for (int x = 0; x < local_nx; x++)
                    {
                        int global_idx = (offset_z + z) * ny * nx + (offset_y + y) * nx + (offset_x + x);
                        int local_idx = z * local_ny * local_nx + y * local_nx + x;
                        local_data[local_idx * nc + t] = time_step_data[global_idx];
                    }
                }
            }

            // Send to other ranks
            for (int dest = 1; dest < size; dest++)
            {
                int dest_coord_x = dest % px;
                int dest_coord_y = (dest / px) % py;
                int dest_coord_z = dest / (px * py);

                int dest_local_nx = nx / px + ((dest_coord_x < nx % px) ? 1 : 0);
                int dest_local_ny = ny / py + ((dest_coord_y < ny % py) ? 1 : 0);
                int dest_local_nz = nz / pz + ((dest_coord_z < nz % pz) ? 1 : 0);

                int dest_offset_x = (nx / px) * dest_coord_x + ((dest_coord_x < nx % px) ? dest_coord_x : nx % px);
                int dest_offset_y = (ny / py) * dest_coord_y + ((dest_coord_y < ny % py) ? dest_coord_y : ny % py);
                int dest_offset_z = (nz / pz) * dest_coord_z + ((dest_coord_z < nz % pz) ? dest_coord_z : nz % pz);

                float *dest_buffer = (float *)malloc(dest_local_nx * dest_local_ny * dest_local_nz * sizeof(float));

                for (int z = 0; z < dest_local_nz; z++)
                {
                    for (int y = 0; y < dest_local_ny; y++)
                    {
                        for (int x = 0; x < dest_local_nx; x++)
                        {
                            int global_idx = (dest_offset_z + z) * ny * nx + (dest_offset_y + y) * nx + (dest_offset_x + x);
                            int buf_idx = z * dest_local_ny * dest_local_nx + y * dest_local_nx + x;
                            dest_buffer[buf_idx] = time_step_data[global_idx];
                        }
                    }
                }

                MPI_Send(dest_buffer, dest_local_nx * dest_local_ny * dest_local_nz, MPI_FLOAT, dest, t, MPI_COMM_WORLD);
                free(dest_buffer);
            }

            free(time_step_data);
        }
        else
        {
            // Receive data from rank 0
            MPI_Recv(&local_data[t], local_nx * local_ny * local_nz, MPI_FLOAT, 0, t, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // Free global data after distribution
    if (rank == 0 && global_data)
    {
        free(global_data);
    }

    // End of file read and data distribution timing
    time2 = MPI_Wtime();
    read_time = time2 - time1;
    MPI_Reduce(&read_time, &max_read_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Process and exchange ghost cells

    // Allocate memory for ghost layers
    // 6 faces: -x, +x, -y, +y, -z, +z
    // Calculate sizes for each face
    int ghost_sizes[6];
    ghost_sizes[0] = local_ny * local_nz; // -x face
    ghost_sizes[1] = local_ny * local_nz; // +x face
    ghost_sizes[2] = local_nx * local_nz; // -y face
    ghost_sizes[3] = local_nx * local_nz; // +y face
    ghost_sizes[4] = local_nx * local_ny; // -z face
    ghost_sizes[5] = local_nx * local_ny; // +z face

    int total_ghost_size = 0;
    for (int i = 0; i < 6; i++)
    {
        total_ghost_size += ghost_sizes[i];
    }

    float *ghost_layers = (float *)malloc(total_ghost_size * nc * sizeof(float));
    if (!ghost_layers)
    {
        fprintf(stderr, "Error: Memory allocation failed for ghost layers on rank %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // For each time step
    for (int t = 0; t < nc; t++)
    {
        // Calculate neighbor ranks
        int neighbor_ranks[6];
        neighbor_ranks[0] = (coord_x > 0) ? rank - 1 : MPI_PROC_NULL;            // -x neighbor
        neighbor_ranks[1] = (coord_x < px - 1) ? rank + 1 : MPI_PROC_NULL;       // +x neighbor
        neighbor_ranks[2] = (coord_y > 0) ? rank - px : MPI_PROC_NULL;           // -y neighbor
        neighbor_ranks[3] = (coord_y < py - 1) ? rank + px : MPI_PROC_NULL;      // +y neighbor
        neighbor_ranks[4] = (coord_z > 0) ? rank - px * py : MPI_PROC_NULL;      // -z neighbor
        neighbor_ranks[5] = (coord_z < pz - 1) ? rank + px * py : MPI_PROC_NULL; // +z neighbor

        // Extract and send boundary data to neighboring processes
        // For each face direction
        MPI_Request send_requests[6], recv_requests[6];

        for (int dir = 0; dir < 6; dir++)
        {
            if (neighbor_ranks[dir] != MPI_PROC_NULL)
            {
                // Allocate buffer for sending
                float *send_buffer = (float *)malloc(ghost_sizes[dir] * sizeof(float));

                // Extract boundary data
                if (dir == 0)
                { // -x face, send leftmost plane
                    for (int z = 0; z < local_nz; z++)
                    {
                        for (int y = 0; y < local_ny; y++)
                        {
                            send_buffer[z * local_ny + y] = local_data[(z * local_ny * local_nx + y * local_nx + 0) * nc + t];
                        }
                    }
                }
                else if (dir == 1)
                { // +x face, send rightmost plane
                    for (int z = 0; z < local_nz; z++)
                    {
                        for (int y = 0; y < local_ny; y++)
                        {
                            send_buffer[z * local_ny + y] = local_data[(z * local_ny * local_nx + y * local_nx + (local_nx - 1)) * nc + t];
                        }
                    }
                }
                else if (dir == 2)
                { // -y face, send bottom plane
                    for (int z = 0; z < local_nz; z++)
                    {
                        for (int x = 0; x < local_nx; x++)
                        {
                            send_buffer[z * local_nx + x] = local_data[(z * local_ny * local_nx + 0 * local_nx + x) * nc + t];
                        }
                    }
                }
                else if (dir == 3)
                { // +y face, send top plane
                    for (int z = 0; z < local_nz; z++)
                    {
                        for (int x = 0; x < local_nx; x++)
                        {
                            send_buffer[z * local_nx + x] = local_data[(z * local_ny * local_nx + (local_ny - 1) * local_nx + x) * nc + t];
                        }
                    }
                }
                else if (dir == 4)
                { // -z face, send front plane
                    for (int y = 0; y < local_ny; y++)
                    {
                        for (int x = 0; x < local_nx; x++)
                        {
                            send_buffer[y * local_nx + x] = local_data[(0 * local_ny * local_nx + y * local_nx + x) * nc + t];
                        }
                    }
                }
                else if (dir == 5)
                { // +z face, send back plane
                    for (int y = 0; y < local_ny; y++)
                    {
                        for (int x = 0; x < local_nx; x++)
                        {
                            send_buffer[y * local_nx + x] = local_data[((local_nz - 1) * local_ny * local_nx + y * local_nx + x) * nc + t];
                        }
                    }
                }

                // Calculate offset in ghost_layers for receiving
                int ghost_offset = 0;
                for (int i = 0; i < dir; i++)
                {
                    ghost_offset += ghost_sizes[i];
                }

                // Send and receive
                MPI_Isend(send_buffer, ghost_sizes[dir], MPI_FLOAT, neighbor_ranks[dir], 100 + dir, MPI_COMM_WORLD, &send_requests[dir]);
                MPI_Irecv(&ghost_layers[(ghost_offset)*nc + t], ghost_sizes[dir], MPI_FLOAT, neighbor_ranks[dir],
                          100 + ((dir % 2 == 0) ? dir + 1 : dir - 1), MPI_COMM_WORLD, &recv_requests[dir]);

                free(send_buffer);
            }
            else
            {
                send_requests[dir] = MPI_REQUEST_NULL;
                recv_requests[dir] = MPI_REQUEST_NULL;
            }
        }

        // Wait for all communications to complete
        MPI_Waitall(6, send_requests, MPI_STATUSES_IGNORE);
        MPI_Waitall(6, recv_requests, MPI_STATUSES_IGNORE);
    }

    // Process the local domain for each time step
    int *local_min_count = (int *)calloc(nc, sizeof(int));
    int *local_max_count = (int *)calloc(nc, sizeof(int));
    float *local_min_value = (float *)malloc(nc * sizeof(float));
    float *local_max_value = (float *)malloc(nc * sizeof(float));

    // Initialize min/max values
    for (int t = 0; t < nc; t++)
    {
        local_min_value[t] = FLT_MAX;
        local_max_value[t] = -FLT_MAX;
    }

    // Process all points in the local domain
    for (int z = 0; z < local_nz; z++)
    {
        for (int y = 0; y < local_ny; y++)
        {
            for (int x = 0; x < local_nx; x++)
            {
                for (int t = 0; t < nc; t++)
                {
                    float value = local_data[(z * local_ny * local_nx + y * local_nx + x) * nc + t];

                    // Update global min/max
                    if (value < local_min_value[t])
                    {
                        local_min_value[t] = value;
                    }
                    if (value > local_max_value[t])
                    {
                        local_max_value[t] = value;
                    }

                    // Check for local minima and maxima
                    if (isLocalMinimum(local_data, x, y, z, t, nx, ny, nz, nc, local_nx, local_ny, local_nz,
                                       offset_x, offset_y, offset_z, ghost_layers))
                    {
                        local_min_count[t]++;
                    }

                    if (isLocalMaximum(local_data, x, y, z, t, nx, ny, nz, nc, local_nx, local_ny, local_nz,
                                       offset_x, offset_y, offset_z, ghost_layers))
                    {
                        local_max_count[t]++;
                    }
                }
            }
        }
    }

    // Gather results from all processes
    int *global_min_count = NULL;
    int *global_max_count = NULL;
    float *global_min_value = NULL;
    float *global_max_value = NULL;

    if (rank == 0)
    {
        global_min_count = (int *)malloc(nc * sizeof(int));
        global_max_count = (int *)malloc(nc * sizeof(int));
        global_min_value = (float *)malloc(nc * sizeof(float));
        global_max_value = (float *)malloc(nc * sizeof(float));
    }

    // Reduce min/max counts
    MPI_Reduce(local_min_count, global_min_count, nc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_max_count, global_max_count, nc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reduce min/max values
    MPI_Reduce(local_min_value, global_min_value, nc, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_max_value, global_max_value, nc, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    // End of main code timing
    time3 = MPI_Wtime();
    main_time = time3 - time2;
    MPI_Reduce(&main_time, &max_main_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Write results to output file
    if (rank == 0)
    {
        FILE *outfile = fopen(output_file, "w");
        if (!outfile)
        {
            fprintf(stderr, "Error opening output file: %s\n", output_file);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Line 1: (count of local minima, count of local maxima), ...
        for (int t = 0; t < nc; t++)
        {
            fprintf(outfile, "(%d, %d)", global_min_count[t], global_max_count[t]);
            if (t < nc - 1)
            {
                fprintf(outfile, ", ");
            }
        }
        fprintf(outfile, "\n");

        // Line 2: (global minimum, global maximum), ...
        for (int t = 0; t < nc; t++)
        {
            fprintf(outfile, "(%f, %f)", global_min_value[t], global_max_value[t]);
            if (t < nc - 1)
            {
                fprintf(outfile, ", ");
            }
        }
        fprintf(outfile, "\n");

        // End of total time
        time4 = MPI_Wtime();
        total_time = time4 - time1;
        MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // Print timing information
        fprintf(outfile, "Read time: %f seconds\n", max_read_time);
        fprintf(outfile, "Main processing time: %f seconds\n", max_main_time);
        fprintf(outfile, "Total execution time: %f seconds\n", max_total_time);

        fclose(outfile);

        // Free allocated memory for global results
        free(global_min_count);
        free(global_max_count);
        free(global_min_value);
        free(global_max_value);
    }

    // Free locally allocated memory
    free(local_data);
    free(ghost_layers);
    free(local_min_count);
    free(local_max_count);
    free(local_min_value);
    free(local_max_value);

    // Clean up MPI datatype
    MPI_Type_free(&local_domain_type);

    MPI_Finalize();
    return 0;
}