#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.hpp"
#include <string.h>

int nx, ny;
double *C;
double *r, *p, *Ap;
double dt, dx, dy;
double diag_coeff;
double offdiag_coeff;
double D;

void compute_matrix_vector_product(double *x, double *result)
{
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            int idx = i * nx + j;
            result[idx] = diag_coeff * x[idx];

            // Interior points use central differences
            if (j > 0)
                result[idx] += offdiag_coeff * x[i * nx + (j - 1)];
            if (j < nx - 1)
                result[idx] += offdiag_coeff * x[i * nx + (j + 1)];
            if (i > 0)
                result[idx] += offdiag_coeff * x[(i - 1) * nx + j];
            if (i < ny - 1)
                result[idx] += offdiag_coeff * x[(i + 1) * nx + j];

            // Boundaries use the ghost cell values (zero-flux)
            if (j == 0)
                result[idx] += offdiag_coeff * x[idx];
            if (j == nx - 1)
                result[idx] += offdiag_coeff * x[idx];
            if (i == 0)
                result[idx] += offdiag_coeff * x[idx];
            if (i == ny - 1)
                result[idx] += offdiag_coeff * x[idx];
        }
    }
}

void init(double *C0, double length, double width, int nx_, int ny_,
          double dt_, double D, int rank, int num_procs)
{
    nx = nx_;
    ny = ny_;
    C = C0;
    dx = length / (nx - 1.0);
    dy = width / (ny - 1.0);
    dt = dt_;

    C = (double *)calloc(nx * ny, sizeof(double));
    r = (double *)calloc(nx * ny, sizeof(double));
    p = (double *)calloc(nx * ny, sizeof(double));
    Ap = (double *)calloc(nx * ny, sizeof(double));

    memcpy(C, C0, nx * ny * sizeof(double));

    diag_coeff = 1.0 + 4.0 * D * dt / (dx * dx);
    offdiag_coeff = -D * dt / (dx * dx);
}

void solve_linear_system()
{
    const double tolerance = 1e-7;
    const int max_iter = 1000;

    compute_matrix_vector_product(C, Ap);
    double r_norm = 0.0;
    for (int i = 0; i < nx * ny; i++)
    {
        r[i] = C[i] - Ap[i];
        p[i] = r[i];
        r_norm += r[i] * r[i];
    }

    int iter = 0;
    while (iter < max_iter && r_norm > tolerance * tolerance)
    {
        compute_matrix_vector_product(p, Ap);

        double pAp = 0.0;
        for (int i = 0; i < nx * ny; i++)
        {
            pAp += p[i] * Ap[i];
        }
        double alpha = r_norm / pAp;

        double r_norm_new = 0.0;
        for (int i = 0; i < nx * ny; i++)
        {
            C[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            r_norm_new += r[i] * r[i];
        }

        double beta = r_norm_new / r_norm;
        for (int i = 0; i < nx * ny; i++)
        {
            p[i] = r[i] + beta * p[i];
        }

        r_norm = r_norm_new;
        iter++;
    }
}

void step()
{
    solve_linear_system();
}

void transfer(double *result)
{
    if (result != C)
    {
        memcpy(result, C, nx * ny * sizeof(double));
    }
}

void free_memory()
{
    free(C);
    free(r);
    free(p);
    free(Ap);
}