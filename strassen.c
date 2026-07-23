#include <sys/time.h>
#include "strassen.h"

#include <stdlib.h>


extern unsigned short int CUTOFF;


double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double) tv.tv_usec / 1000000;
    return seconds;
}


void populate_matrix(matrix A) {
    srand(10);
    double rand_max = 1;
    double rand_min = 0;
    double rand_range = rand_max - rand_min;
    for (int i = 0; i < A.N; i++) {
        for (int j = 0; j < A.N; j++) {
            matrix(A, i, j) = rand_min + rand_range * ((double) rand() / RAND_MAX);
        }
    }
}

// Test if the add and subtract functions are vectorized differently
void matrix_add(const matrix *restrict A, const matrix *restrict B, matrix *restrict C) {
    const int N = A->N;
    // Perform element-wise addition for entire matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix((*C), i, j) = matrix((*A), i, j) + matrix((*B), i, j);
        }
    }
}


void matrix_subtract(const matrix *restrict A, const matrix *restrict B, matrix *restrict C) {
    const int N = A->N;
    // Perform element-wise subtraction for entire matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix((*C), i, j) = matrix((*A), i, j) - matrix((*B), i, j);
        }
    }
}


void naive_multiply(const matrix *restrict A, const matrix *restrict B, matrix *restrict C) {
    const int N = A->N;
    // Naive matrix multiplication in the order i-k-j for cache locality
    // *restrict keyword also help the optimizer assume to overlapping memory addresses
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) matrix((*C), i, j) = 0.0;

        for (int k = 0; k < N; k++) {
            const double a = matrix((*A), i, k);  // Lifted out of the inner loop
            for (int j = 0; j < N; j++) {
                matrix((*C), i, j) += a * matrix((*B), k, j);
            }
        }
    }
}


void strassen_multiply(const matrix *restrict A, const matrix *restrict B, matrix *restrict C) {
    // Use naive matrix multiplication if the CUTOFF is reached
    if (A->N <= CUTOFF) {
        naive_multiply(A, B, C);
        return;
    }

    int N = A->N;
    int mid = N / 2;

    // Partition the matrices into submatrix views
    matrix A11 = { .data = A->data, .N = mid, .tda = A->tda };
    matrix A12 = { .data = A->data + mid, .N = mid, .tda = A->tda };
    matrix A21 = { .data = A->data + mid * A->tda, .N = mid, .tda = A->tda };
    matrix A22 = { .data = A->data + mid * A->tda + mid, .N = mid, .tda = A->tda };

    matrix B11 = { .data = B->data, .N = mid, .tda = B->tda };
    matrix B12 = { .data = B->data + mid, .N = mid, .tda = B->tda };
    matrix B21 = { .data = B->data + mid * B->tda, .N = mid, .tda = B->tda };
    matrix B22 = { .data = B->data + mid * B->tda + mid, .N = mid, .tda = B->tda };

    matrix C11 = { .data = C->data, .N = mid, .tda = C->tda };
    matrix C12 = { .data = C->data + mid, .N = mid, .tda = C->tda };
    matrix C21 = { .data = C->data + mid * C->tda, .N = mid, .tda = C->tda };
    matrix C22 = { .data = C->data + mid * C->tda + mid, .N = mid, .tda = C->tda };


}