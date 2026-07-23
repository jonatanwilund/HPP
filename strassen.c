#include <sys/time.h>
#include "strassen.h"
#include <stdlib.h>

extern unsigned short int CUTOFF;
extern unsigned short int MAX_RECURSION_DEPTH;

double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    const double seconds = tv.tv_sec + (double) tv.tv_usec / 1000000;
    return seconds;
}

void populate_matrix(matrix A) {
    srand(10);
    constexpr double rand_max = 1;
    constexpr double rand_min = 0;
    constexpr double rand_range = rand_max - rand_min;
    for (int i = 0; i < A.N; i++) {
        for (int j = 0; j < A.N; j++) {
            matrix(A, i, j) = rand_min + rand_range * ((double) rand() / RAND_MAX);
        }
    }
}

// Test if the add and subtract functions are vectorized differently
void matrix_add(matrix C, matrix A, matrix B) {
    const unsigned int N = A.N;
    double *restrict c_data = C.data;
    const double *restrict a_data = A.data;
    const double *restrict b_data = B.data;
    // Perform element-wise addition for entire matrix
    for (unsigned int i = 0; i < N; i++) {
        double *restrict c_row = c_data + i * C.tda;
        const double *restrict a_row = a_data + i * A.tda;
        const double *restrict b_row = b_data + i * B.tda;
        for (unsigned int j = 0; j < N; j++) {
            c_row[j] = a_row[j] + b_row[j];
        }
    }
}

void matrix_subtract(matrix C, matrix A, matrix B) {
    const unsigned int N = A.N;
    double *restrict c_data = C.data;
    const double *restrict a_data = A.data;
    const double *restrict b_data = B.data;
    // Perform element-wise subtraction for entire matrix
    for (unsigned int i = 0; i < N; i++) {
        double *restrict c_row = c_data + i * C.tda;
        const double *restrict a_row = a_data + i * A.tda;
        const double *restrict b_row = b_data + i * B.tda;
        for (unsigned int j = 0; j < N; j++) {
            c_row[j] = a_row[j] - b_row[j];
        }
    }
}

void assemble_matrix(
    matrix C,
    matrix C11,
    matrix C12,
    matrix C21,
    matrix C22) {
    const unsigned int mid = C11.N;
    for (unsigned int i = 0; i < mid; i++) {
        for (unsigned int j = 0; j < mid; j++) {
            matrix(C, i, j) = matrix(C11, i, j);
            matrix(C, i, j + mid) = matrix(C12, i, j);
            matrix(C, i + mid, j) = matrix(C21, i, j);
            matrix(C, i + mid, j + mid) = matrix(C22, i, j);
        }
    }
}

void naive_multiply(matrix C, matrix A, matrix B) {
    const unsigned int N = A.N;
    double *restrict c_data = C.data;
    const double *restrict a_data = A.data;
    const double *restrict b_data = B.data;
    // Naive matrix multiplication in the order i-k-j for cache locality
    // *restrict keyword also help the optimizer assume to overlapping memory addresses
    for (unsigned int i = 0; i < N; i++) {
        double *restrict c_row = c_data + i * C.tda;
        for (unsigned int j = 0; j < N; j++) c_row[j] = 0.0;

        for (unsigned int k = 0; k < N; k++) {
            const double a = a_data[i * A.tda + k];  // Lifted out of the inner loop
            const double *restrict b_row = b_data + k * B.tda;
            for (unsigned int j = 0; j < N; j++) {
                c_row[j] += a * b_row[j];
            }
        }
    }
}

void strassen_multiply(matrix C, matrix A, matrix B) {
    // Use naive matrix multiplication if the CUTOFF is reached
    if (A.N <= CUTOFF) {
        naive_multiply(C, A, B);
        return;
    }

    unsigned int N = A.N;
    unsigned int mid = N / 2;
    unsigned int num_sub_elements = mid * mid;

    // Partition the matrices into submatrix views
    const matrix A11 = { .data = A.data, .N = mid, .tda = A.tda };
    const matrix A12 = { .data = A.data + mid, .N = mid, .tda = A.tda };
    const matrix A21 = { .data = A.data + mid * A.tda, .N = mid, .tda = A.tda };
    const matrix A22 = { .data = A.data + mid * A.tda + mid, .N = mid, .tda = A.tda };

    const matrix B11 = { .data = B.data, .N = mid, .tda = B.tda };
    const matrix B12 = { .data = B.data + mid, .N = mid, .tda = B.tda };
    const matrix B21 = { .data = B.data + mid * B.tda, .N = mid, .tda = B.tda };
    const matrix B22 = { .data = B.data + mid * B.tda + mid, .N = mid, .tda = B.tda };

    matrix C11 = { .data = C.data, .N = mid, .tda = C.tda };
    matrix C12 = { .data = C.data + mid, .N = mid, .tda = C.tda };
    matrix C21 = { .data = C.data + mid * C.tda, .N = mid, .tda = C.tda };
    matrix C22 = { .data = C.data + mid * C.tda + mid, .N = mid, .tda = C.tda };

    // Allocate one contiguous block of memory for all sub-problems
    double *mem_block = malloc(17 * num_sub_elements * sizeof(double));  // 7 mult. results, 10 add./subtr. results

    // Partition the memory block into views
    matrix M[7];
    matrix S[5];
    matrix T[5];
    for (int i = 0; i < 7; i++) {
        M[i] = (matrix){ .data = mem_block + i * num_sub_elements, .N = mid, .tda = mid };
    }
    for (int i = 0; i < 5; i++) {
        S[i] = (matrix){ .data = mem_block + (7 + i) * num_sub_elements, .N = mid, .tda = mid };
        T[i] = (matrix){ .data = mem_block + (12 + i) * num_sub_elements, .N = mid, .tda = mid };
    }

    // M1 = (A11 + A22) * (B11 + B22)
#pragma omp task depend(out: M[0])
    {
        matrix_add(S[0], A11, A22);
        matrix_add(T[0], B11, B22);
        strassen_multiply(M[0], S[0], T[0]);
    }

    // M2 = (A21 + A22) * B11    Note: T2 = B11, no T buffer needed
#pragma omp task depend(out: M[1])
    {
        matrix_add(S[1], A21, A22);
        strassen_multiply(M[1], S[1], B11);
    }

    // M3 = A11 * (B12 - B22)    Note: S2 = A11, no S buffer needed
#pragma omp task depend(out: M[2])
    {
        matrix_subtract(T[1], B12, B22);
        strassen_multiply(M[2], A11, T[1]);
    }

    // M4 = A22 * (B21 - B11)    Note: S4 = A22 directly, no S buffer needed
#pragma omp task depend(out: M[3])
    {
        matrix_subtract(T[2], B21, B11);
        strassen_multiply(M[3], A22, T[2]);
    }

    // M5 = (A11 + A12) * B22    Note: T5 = B22 directly, no T buffer needed
#pragma omp task depend(out: M[4])
    {
        matrix_add(S[2], A11, A12);
        strassen_multiply(M[4], S[2], B22);
    }

    // M6 = (A21 - A11) * (B11 + B12)
#pragma omp task depend(out: M[5])
    {
        matrix_subtract(S[3], A21, A11);
        matrix_add(T[3], B11, B12);
        strassen_multiply(M[5], S[3], T[3]);
    }

    // M7 = (A12 - A22) * (B21 + B22)
#pragma omp task depend(out: M[6])
    {
        matrix_subtract(S[4], A12, A22);
        matrix_add(T[4], B21, B22);
        strassen_multiply(M[6], S[4], T[4]);
    }

    // Assemble temporary matrix products in matrix C's sub-matrices
#pragma omp task depend(in: M[0], M[3], M[4], M[6]) depend(out: C11)
    {
        double *restrict c11_data = C11.data;
        const double *restrict m1_data = M[0].data;
        const double *restrict m4_data = M[3].data;
        const double *restrict m5_data = M[4].data;
        const double *restrict m7_data = M[6].data;
        for (int i = 0; i < mid; i++) {
            double *restrict c11_row = c11_data + i * C11.tda;
            const double *restrict m1_row = m1_data + i * M[0].tda;
            const double *restrict m4_row = m4_data + i * M[3].tda;
            const double *restrict m5_row = m5_data + i * M[4].tda;
            const double *restrict m7_row = m7_data + i * M[6].tda;
            for (int j = 0; j < mid; j++) {
                c11_row[j] = m1_row[j] + m4_row[j] - m5_row[j] + m7_row[j];
            }
        }
    }

#pragma omp task depend(in: M[2], M[4]) depend(out: C12)
    matrix_add(C12, M[2], M[4]);

#pragma omp task depend(in: M[1], M[3]) depend(out: C21)
    matrix_add(C21, M[1], M[3]);

#pragma omp task depend(in: M[0], M[1], M[2], M[5]) depend(out: C22)
    {
        double *restrict c22_data = C22.data;
        const double *restrict m1_data = M[0].data;
        const double *restrict m2_data = M[1].data;
        const double *restrict m3_data = M[2].data;
        const double *restrict m6_data = M[5].data;
        for (int i = 0; i < mid; i++) {
            double *restrict c22_row = c22_data + i * C22.tda;
            const double *restrict m1_row = m1_data + i * M[0].tda;
            const double *restrict m2_row = m2_data + i * M[1].tda;
            const double *restrict m3_row = m3_data + i * M[2].tda;
            const double *restrict m6_row = m6_data + i * M[5].tda;
            for (int j = 0; j < mid; j++) {
                c22_row[j] = m1_row[j] - m2_row[j] + m3_row[j] + m6_row[j];
            }
        }
    }

#pragma omp taskwait

    free(mem_block);
}