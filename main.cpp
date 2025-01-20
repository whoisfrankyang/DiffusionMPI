#include "solver.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <chrono>
#include <cstring>

#ifdef MPI
#include <mpi.h>
#endif

// void integrate_output_files(const char *output_folder, const char *integrated_file, int num_steps, int nx, int ny)
// {
//     char filename[200];
//     FILE *out_fp = fopen(integrated_file, "w");
//     if (out_fp == NULL)
//     {
//         printf("Error: Could not open integrated output file %s\n", integrated_file);
//         return;
//     }

//     fprintf(out_fp, "# Integrated simulation output\n");
//     fprintf(out_fp, "# Grid dimensions: %d x %d\n", nx, ny);
//     fprintf(out_fp, "# Number of timesteps: %d\n\n", num_steps);

//     // Integrate initial state
//     sprintf(filename, "%s/initial_state.txt", output_folder);
//     FILE *fp = fopen(filename, "r");
//     if (fp != NULL)
//     {
//         fprintf(out_fp, "# Initial state\n");
//         char line[1024];
//         while (fgets(line, sizeof(line), fp))
//         {
//             fprintf(out_fp, "%s", line);
//         }
//         fclose(fp);
//         fprintf(out_fp, "\n");
//     }

//     // Integrate each timestep
//     for (int timestep = 0; timestep < num_steps; timestep++)
//     {
//         sprintf(filename, "%s/state_%03d.txt", output_folder, timestep);
//         fp = fopen(filename, "r");
//         if (fp != NULL)
//         {
//             fprintf(out_fp, "# Timestep %d\n", timestep);
//             char line[1024];
//             while (fgets(line, sizeof(line), fp))
//             {
//                 fprintf(out_fp, "%s", line);
//             }
//             fclose(fp);
//             fprintf(out_fp, "\n");
//         }
//     }

//     // Integrate final state
//     sprintf(filename, "%s/final_state.txt", output_folder);
//     fp = fopen(filename, "r");
//     if (fp != NULL)
//     {
//         fprintf(out_fp, "# Final state\n");
//         char line[1024];
//         while (fgets(line, sizeof(line), fp))
//         {
//             fprintf(out_fp, "%s", line);
//         }
//         fclose(fp);
//     }

//     fclose(out_fp);
//     printf("Integrated output written to %s\n", integrated_file);
// }

// void save_to_file(const char *filename, double *data, int nx, int ny)
// {
//     FILE *fp = fopen(filename, "w");
//     if (fp == NULL)
//     {
//         printf("Error: Could not open file %s\n", filename);
//         return;
//     }

//     double max_val = data[0];
//     double min_val = data[0];

//     for (int i = 0; i < nx; i++)
//     {
//         for (int j = 0; j < ny; j++)
//         {
//             fprintf(fp, "%.6f ", data[i * ny + j]);
//             if (data[i * ny + j] > max_val)
//                 max_val = data[i * ny + j];
//             if (data[i * ny + j] < min_val)
//                 min_val = data[i * ny + j];
//         }
//         fprintf(fp, "\n");
//     }
//     fclose(fp);

//     printf("File %s written. Value range: [%f, %f]\n", filename, min_val, max_val);
// }

int main(int argc, char **argv)
{
    int nx = 50;
    int ny = 50;
    double length = 1.0;
    double width = 1.0;
    double dt = 0.001;
    int num_iterations = 200;
    double D = 0.1;

    int rank = 0,
        num_procs = 1;

#ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    int cur_arg = 1;
    int num_args = argc - 1;
    char output_file[256];
    bool output = false;

    while (num_args > 0)
    {
        if (num_args == 1)
        {
            fprintf(stderr, "Missing argument value for %s\n", argv[cur_arg]);
            return 1;
        }

        if (strcmp(argv[cur_arg], "--nx") == 0)
        {
            nx = atoi(argv[cur_arg + 1]);
        }
        else if (strcmp(argv[cur_arg], "--ny") == 0)
        {
            ny = atoi(argv[cur_arg + 1]);
        }
        else if (strcmp(argv[cur_arg], "--iter") == 0)
        {
            num_iterations = atoi(argv[cur_arg + 1]);
        }
        else if (strcmp(argv[cur_arg], "--output") == 0)
        {
            strcpy(output_file, argv[cur_arg + 1]);
            output = true;
        }
        else
        {
            fprintf(stderr, "Unknown argument: %s\n", argv[cur_arg]);
            return 1;
        }

        cur_arg += 2;
        num_args -= 2;
    }

    FILE *fptr;

    if (rank == 0 && output)
    {
        fptr = fopen(output_file, "w");
        fwrite(&nx, sizeof(int), 1, fptr);
        fwrite(&ny, sizeof(int), 1, fptr);
        fwrite(&num_iterations, sizeof(int), 1, fptr);
    }

    double cx = length / 2;
    double cy = width / 2;
    double sigma = 0.1;

    double *C0 = nullptr;
    if (rank == 0)
    {
        C0 = (double *)calloc(nx * ny, sizeof(double));
        double peak_scale = 10.0; // Scale factor to make the peak higher

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                double x = length * i / (nx - 1.0);
                double y = width * j / (ny - 1.0);
                double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                C0[i * ny + j] = peak_scale * exp(-r2 / (2 * sigma * sigma)); // Amplify peak
            }
        }
    }

    clock_t init_start = clock();
    init(C0, length, width, nx, ny, dt, D, rank, num_procs);
    clock_t init_end = clock();
    fprintf(stderr, "Initialization time for rank %d: %f\n", rank, (double)(init_end - init_start) / CLOCKS_PER_SEC);

    clock_t start = clock();
    // printf("Starting simulation...\n");
    int timestep = 0;
    for (; timestep < num_iterations; timestep++)
    {
        if (output)
        {
            transfer(C0);
            if (rank == 0)
            {
                fwrite(C0, sizeof(double), nx * ny, fptr);
            }
        }
        step();
    }

    if (output)
    {
        transfer(C0);
        if (rank == 0)
        {
            fwrite(C0, sizeof(double), nx * ny, fptr);
        }
    }

    transfer(C0);
    clock_t end = clock();
    fprintf(stderr, "Execution time for rank %d: %f\n", rank, (double)(end - start) / CLOCKS_PER_SEC);

    free_memory();
    if (rank == 0)
    {
        free(C0);
        if (output)
        {
            fclose(fptr);
        }
    }

#ifdef MPI
    MPI_Finalize();
#endif
    return 0;
}
