#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "strassen.h"

unsigned short int CUTOFF = 1;
unsigned short int MAX_RECURSION_DEPTH = 4;

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("Usage: %s <N> <N_cutoff> <max_recursion_depth> <n_threads>\n", argv[0]);
        return -1;
    }
    const unsigned int N = atoi(argv[1]);  // Matrix size
    CUTOFF = atoi(argv[2]);  // Size when the recursion should stop
    MAX_RECURSION_DEPTH = atoi(argv[3]);  // Max number of recursion levels before stopping creating omp tasks
    const unsigned int n_threads = atoi(argv[4]);

    if (N == 0 || (N & (N - 1)) != 0) {
        printf("Matrix size must be a power of two.\n");
        return -1;
    }

    omp_set_num_threads(n_threads);

    // Allocate memory for A, B, and C
    double *A_data = malloc((N * N) * sizeof(double));
    double *B_data = malloc((N * N) * sizeof(double));
    double *C_data = malloc((N * N) * sizeof(double));

    // Construct matrices to multiply
    matrix A = { .data = A_data, .N = N, .tda = N };
    matrix B = { .data = B_data, .N = N, .tda = N };
    matrix C = { .data = C_data, .N = N, .tda = N };

    populate_matrix(A);
    populate_matrix(B);

    double start_time = get_wall_seconds();

#pragma omp parallel
    {
#pragma omp single
        {
            strassen_multiply(C, A, B);
        }
    }

    double elapsed_time = get_wall_seconds() - start_time;

    printf("N=%u, Threads=%u, Cutoff=%u, Time=%.4f s\n", N, n_threads, CUTOFF, elapsed_time);

    // Cleanup
    free(A_data);
    free(B_data);
    free(C_data);
    return 0;
}