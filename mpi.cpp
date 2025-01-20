#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include "solver.hpp"

#define C(i, j) C[(i) * ny + (j)]
#define result(i, j) result[(i) * ny + (j)]
#define x(i, j) x[(i) * ny + (j)]
#define p(i, j) p[(i) * ny + (j)]
#define r(i, j) r[(i) * ny + (j)]
#define Ap(i, j) Ap[(i) * ny + (j)]

int nx, ny;
double dx, dy, dt;
double *C, *C_new, *r, *p, *Ap;
double diag_coeff, offdiag_coeff;
int *exchangeCounts, *displs;
double D;
int rank, num_procs, local_nx;

void init(double *C0, double length, double width, int nx_, int ny_, double dt_, double D, int _rank, int _num_procs)
{
    nx = nx_;
    ny = ny_;
    dt = dt_;

    dx = length / (nx - 1.0);
    dy = width / (ny - 1.0);

    rank = _rank;
    num_procs = _num_procs;

    local_nx = nx / num_procs;
    if (rank == num_procs - 1)
    {
        local_nx += nx % num_procs;
    }
    int local_size = (local_nx + 2) * ny;

    C = (double *)calloc(local_size, sizeof(double));
    r = (double *)calloc(local_size, sizeof(double));
    p = (double *)calloc(local_size, sizeof(double));
    Ap = (double *)calloc(local_size, sizeof(double));

    exchangeCounts = (int *)calloc(num_procs, sizeof(int));
    displs = (int *)calloc(num_procs, sizeof(int));

    for (int i = 0; i < num_procs; i++)
    {
        exchangeCounts[i] = nx / num_procs;
        if (i == num_procs - 1)
        {
            exchangeCounts[i] += nx % num_procs;
        }
        exchangeCounts[i] *= ny;
        displs[i] = i * (nx / num_procs) * ny;
    }

    MPI_Scatterv(C0, exchangeCounts, displs, MPI_DOUBLE, &C(1, 0), local_nx * ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    diag_coeff = 1.0 + 4.0 * D * dt / (dx * dx);
    offdiag_coeff = -D * dt / (dx * dx);
}

void exchangeCells(double *x)
{
    // Compute ranks of neighbors
    int prev_rank = rank - 1;
    int next_rank = rank + 1;

    // Adjust for boundary conditions
    if (rank == 0)
        prev_rank = MPI_PROC_NULL; // No previous rank for the first process
    if (rank == num_procs - 1)
        next_rank = MPI_PROC_NULL; // No next rank for the last process

    // Exchange with the next rank
    MPI_Sendrecv(&x(local_nx, 0), ny, MPI_DOUBLE, next_rank, 0,     // Send to next
                 &x(local_nx + 1, 0), ny, MPI_DOUBLE, next_rank, 0, // Receive from next
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange with the previous rank
    MPI_Sendrecv(&x(1, 0), ny, MPI_DOUBLE, prev_rank, 0, // Send to previous
                 &x(0, 0), ny, MPI_DOUBLE, prev_rank, 0, // Receive from previous
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void computeMatVecProduct(double *x, double *result)
{
    for (int i = 1; i < local_nx + 1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            result(i, j) = diag_coeff * x(i, j);
            if (j > 0)
                result(i, j) += offdiag_coeff * x(i, j - 1);
            if (j < ny - 1)
                result(i, j) += offdiag_coeff * x(i, j + 1);
            if (rank == 0)
            {
                if (i > 1)
                    result(i, j) += offdiag_coeff * x(i - 1, j);
                if (i <= local_nx)
                    result(i, j) += offdiag_coeff * x(i + 1, j);
            }
            else if (rank == num_procs - 1)
            {
                if (i < local_nx)
                    result(i, j) += offdiag_coeff * x(i + 1, j);
                if (i >= 1)
                    result(i, j) += offdiag_coeff * x(i - 1, j);
            }
            else
            {
                if (i >= 1)
                    result(i, j) += offdiag_coeff * x(i - 1, j);
                if (i <= local_nx)
                    result(i, j) += offdiag_coeff * x(i + 1, j);
            }

            // zero flux
            if (j == 0)
                result(i, j) += offdiag_coeff * x(i, j);
            if (j == ny - 1)
                result(i, j) += offdiag_coeff * x(i, j);
            if (rank == 0 && i == 1)
            {
                result(i, j) += offdiag_coeff * x(i, j);
            }
            else if (rank == num_procs - 1 && i == local_nx)
            {
                result(i, j) += offdiag_coeff * x(i, j);
            }
        }
    }
}

void solveLinearSystem()
{
    const double tolerance = 1e-7;
    const int max_iter = 1000;
    computeMatVecProduct(C, Ap);
    double r_norm = 0.0;
    for (int i = 1; i < local_nx + 1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            r(i, j) = C(i, j) - Ap(i, j);
            p(i, j) = r(i, j);
            r_norm += r(i, j) * r(i, j);
        }
    }

    double all_r_norm;
    MPI_Allreduce(&r_norm, &all_r_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int iter = 0;
    while (iter < max_iter && all_r_norm > tolerance * tolerance)
    {
        exchangeCells(p);
        computeMatVecProduct(p, Ap);

        double pAp = 0.0;
        for (int i = 1; i < local_nx + 1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                pAp += p(i, j) * Ap(i, j);
            }
        }
        double all_pAp;
        MPI_Allreduce(&pAp, &all_pAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double alpha = all_r_norm / all_pAp;

        double r_norm_new = 0.0;
        for (int i = 1; i < local_nx + 1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                C(i, j) += alpha * p(i, j);
                r(i, j) -= alpha * Ap(i, j);
                r_norm_new += r(i, j) * r(i, j);
            }
        }

        double all_r_norm_new;
        MPI_Allreduce(&r_norm_new, &all_r_norm_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double beta = all_r_norm_new / all_r_norm;
        for (int i = 1; i < local_nx + 1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                p(i, j) = r(i, j) + beta * p(i, j);
            }
        }
        all_r_norm = all_r_norm_new;
        iter++;
    }
}

void step()
{
    exchangeCells(C);
    solveLinearSystem();
}

void transfer(double *result)
{
    MPI_Gatherv(&C(1, 0), local_nx * ny, MPI_DOUBLE, result, exchangeCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void free_memory()
{
    free(C);
    free(r);
    free(p);
    free(Ap);
    free(exchangeCounts);
    free(displs);
}