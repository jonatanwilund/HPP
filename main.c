#include <stdio.h>
#include <stdlib.h>
#include "strassen.h"

unsigned short int CUTOFF = 1;

int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Usage: %s <N> <N_cutoff> <n_threads>\n", argv[0]);
        return -1;
    }
    const unsigned int N = atoi(argv[1]);  // Matrix size
    CUTOFF = atoi(argv[2]);  // Size when the recursion should stop
    const unsigned int n_threads = atoi(argv[3]);

    //(x > 0) && ((x & (x - 1)) == 0)
    if ((N > 0) && ((N & (N - 1)) == 0)) {
        printf("Matrix size must be a power of two.\n");
        return -1;
    }

    // Allocate memory for A, B, and C
    double *A = malloc((N * N) * sizeof(double));
    double *B = malloc((N * N) * sizeof(double));
    double *C = malloc((N * N) * sizeof(double));

    // Cleanup
    free(A);
    free(B);
    free(C);
    return 0;
}