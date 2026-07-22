#include <sys/time.h>
#include "strassen.h"

#include <stdlib.h>


double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double) tv.tv_usec / 1000000;
    return seconds;
}


void populate_matrix(double *A, int N) {
    srand(10);
    double rand_max = 1;
    double rand_min = 0;
    double rand_range = rand_max - rand_min;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = rand_min + rand_range * ((double) rand() / RAND_MAX);
        }
    }
}


void matrix_add(
    const double *restrict A,
    const double *restrict B,
    double *restrict C,
    unsigned int N,
    unsigned int stride) {
    // Perform element-wise addition for entire matrix
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
            C[i * stride + j] = A[i * stride + j] + B[i * stride + j];
        }
    }
}


void matrix_subtract(
    const double *restrict A,
    const double *restrict B,
    double *restrict C,
    unsigned int N,
    unsigned int stride) {
    // Perform element-wise subtraction for entire matrix
    const unsigned int NN = N * N;
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
            C[i * stride + j] = A[i * stride + j] - B[i * stride + j];
        }
    }
}


void naive_multiply(const double *restrict A, const double *restrict B, double *restrict C, const unsigned int N) {
    // Naive matrix multiplication in the order i-k-j for cache locality
    // *restrict keyword also help the optimizer assume to overlapping memory addresses
    for (unsigned int i = 0; i < N; i++) {
        double *c = C + i * N;  // Contiguous -> compiler can vectorize well
        for (unsigned int j = 0; j < N; j++) c[j] = 0.0;

        for (unsigned int k = 0; k < N; k++) {
            const double a = A[i * N + k];  // Lifted out of the inner loop
            const double *b = B + k * N; // Contiguous -> compiler can vectorize well
            for (unsigned int j = 0; j < N; j++) {
                c[j] += a * b[j];
            }
        }
    }
}


void strassen_multiply(
    const double *restrict A,
    const double *restrict B,
    double *restrict C,
    unsigned int N,
    unsigned int stride) {
    // Use naive matrix multiplication if the CUTOFF is reached
    if (N <= CUTOFF) {
        naive_multiply(A, B, C, N);
        return;
    }

    unsigned int N_half = N / 2;

    // Partition the matrices
    const double *A11 = A;
    const double *A12 = A + N_half;
    const double *A21 = A + N_half * stride;
    const double *A22 = A + N_half * stride + N_half;
    const double *B11 = B;
    const double *B12 = B + N_half * stride;
    const double *B21 = B + N_half * stride + N_half;
    const double *B22 = B + N_half * stride + N_half;
    double *C11 = C;
    double *C12 = C + N_half * stride;
    double *C21 = C + N_half * stride + N_half;
    double *C22 = C + N_half * stride + N_half;



}