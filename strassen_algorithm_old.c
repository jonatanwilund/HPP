#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}


void populate_matrix(double *A, int N) {
    srand(10);
    double rand_max = 1;
    double rand_min = 0;
    double rand_range = rand_max - rand_min;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = rand_min + rand_range * ((double)rand() / RAND_MAX);
        }
    }
}


// Recursively partition the matrix into two smaller matrices and multiply them
void matmul(double *restrict A, double *restrict B, double *restrict C, const int N, const int stride) {
    if (N < 64) {
        // Initialize the target block of C to 0
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                C[i * stride + j] = 0.0;
            }
        }

        // i-k-j loop order for cache locality
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < N; k++) {
                // Load the A value once per k so we don't look it up repeatedly
                double a_ik = A[i * stride + k];
                for (int j = 0; j < N; j++) {
                    C[i * stride + j] += a_ik * B[k * stride + j];
                }
            }
        }
        return;
    } else {
        int mid = N / 2;

        // One allocation for all 14 sub_matrices that will be multiplied together
        double *sub_matrices_flat = malloc(14 * mid * mid * sizeof(double));
        double *restrict S0 = sub_matrices_flat;
        double *restrict S1 = sub_matrices_flat + mid * mid;
        double *restrict S2 = sub_matrices_flat + 2 * mid * mid;
        double *restrict S3 = sub_matrices_flat + 3 * mid * mid;
        double *restrict S4 = sub_matrices_flat + 4 * mid * mid;
        double *restrict S5 = sub_matrices_flat + 5 * mid * mid;
        double *restrict S6 = sub_matrices_flat + 6 * mid * mid;
        double *restrict S7 = sub_matrices_flat + 7 * mid * mid;
        double *restrict S8 = sub_matrices_flat + 8 * mid * mid;
        double *restrict S9 = sub_matrices_flat + 9 * mid * mid;
        double *restrict S10 = sub_matrices_flat + 10 * mid * mid;
        double *restrict S11 = sub_matrices_flat + 11 * mid * mid;
        double *restrict S12 = sub_matrices_flat + 12 * mid * mid;
        double *restrict S13 = sub_matrices_flat + 13 * mid * mid;

        // One allocation for all 7 "result" matrices
        double *M_flat = malloc(7 * mid * mid * sizeof(double));
        double *restrict M0 = M_flat;
        double *restrict M1 = M_flat + mid * mid;
        double *restrict M2 = M_flat + 2 * mid * mid;
        double *restrict M3 = M_flat + 3 * mid * mid;
        double *restrict M4 = M_flat + 4 * mid * mid;
        double *restrict M5 = M_flat + 5 * mid * mid;
        double *restrict M6 = M_flat + 6 * mid * mid;

        // Partition the matrix into two smaller matrices that will be multiplied
        for (int i = 0; i < mid; i++) {
            // Calculate the starting pointer for the current row of A and B
            double *restrict A11_row = A + (i * stride);
            double *restrict A12_row = A + (i * stride) + mid;
            double *restrict A21_row = A + ((i + mid) * stride);
            double *restrict A22_row = A + ((i + mid) * stride) + mid;

            double *restrict B11_row = B + (i * stride);
            double *restrict B12_row = B + (i * stride) + mid;
            double *restrict B21_row = B + ((i + mid) * stride);
            double *restrict B22_row = B + ((i + mid) * stride) + mid;

            // Calculate the starting pointer for the current row of sub_matrices
            double *restrict S0_row = S0 + (i * mid);
            double *restrict S1_row = S1 + (i * mid);
            double *restrict S2_row = S2 + (i * mid);
            double *restrict S3_row = S3 + (i * mid);
            double *restrict S4_row = S4 + (i * mid);
            double *restrict S5_row = S5 + (i * mid);
            double *restrict S6_row = S6 + (i * mid);
            double *restrict S7_row = S7 + (i * mid);
            double *restrict S8_row = S8 + (i * mid);
            double *restrict S9_row = S9 + (i * mid);
            double *restrict S10_row = S10 + (i * mid);
            double *restrict S11_row = S11 + (i * mid);
            double *restrict S12_row = S12 + (i * mid);
            double *restrict S13_row = S13 + (i * mid);

            // Tell the compiler to ignore potential dependencies
            #pragma GCC ivdep
            for (int j = 0; j < mid; j++) {
                // M1
                S0_row[j] = A11_row[j] + A22_row[j];
                S1_row[j] = B11_row[j] + B22_row[j];

                // M2
                S2_row[j] = A21_row[j] + A22_row[j];
                S3_row[j] = B11_row[j];

                // M3
                S4_row[j] = A11_row[j];
                S5_row[j] = B12_row[j] - B22_row[j];

                // M4
                S6_row[j] = A22_row[j];
                S7_row[j] = B21_row[j] - B11_row[j];

                // M5
                S8_row[j] = A11_row[j] + A12_row[j];
                S9_row[j] = B22_row[j];

                // M6
                S10_row[j] = A21_row[j] - A11_row[j];
                S11_row[j] = B11_row[j] + B12_row[j];

                // M7
                S12_row[j] = A12_row[j] - A22_row[j];
                S13_row[j] = B21_row[j] + B22_row[j];
            }
        }
        // Multiply the submatrices
#pragma omp task shared(S0, S1, M0) firstprivate(mid)
        matmul(S0, S1, M0, mid, mid);

#pragma omp task shared(S2, S3, M1) firstprivate(mid)
        matmul(S2, S3, M1, mid, mid);

#pragma omp task shared(S4, S5, M2) firstprivate(mid)
        matmul(S4, S5, M2, mid, mid);

#pragma omp task shared(S6, S7, M3) firstprivate(mid)
        matmul(S6, S7, M3, mid, mid);

#pragma omp task shared(S8, S9, M4) firstprivate(mid)
        matmul(S8, S9, M4, mid, mid);

#pragma omp task shared(S10, S11, M5) firstprivate(mid)
        matmul(S10, S11, M5, mid, mid);

#pragma omp task shared(S12, S13, M6) firstprivate(mid)
        matmul(S12, S13, M6, mid, mid);

#pragma omp taskwait

        for (int i = 0; i < mid; i++) {
            double *restrict M0_row = M0 + (i * mid);
            double *restrict M1_row = M1 + (i * mid);
            double *restrict M2_row = M2 + (i * mid);
            double *restrict M3_row = M3 + (i * mid);
            double *restrict M4_row = M4 + (i * mid);
            double *restrict M5_row = M5 + (i * mid);
            double *restrict M6_row = M6 + (i * mid);

            double *restrict C11_row = C + (i * stride);
            double *restrict C12_row = C + (i * stride) + mid;
            double *restrict C21_row = C + ((i + mid) * stride);
            double *restrict C22_row = C + ((i + mid) * stride) + mid;
#pragma GCC ivdep
            for (int j = 0; j < mid; j++) {
                C11_row[j] = M0_row[j] + M3_row[j] - M4_row[j] + M6_row[j];  // C11 = M1 + M4 - M5 + M7
                C12_row[j] = M2_row[j] + M4_row[j];  // C12 = M3 + M5
                C21_row[j] = M1_row[j] + M3_row[j];  // C21 = M2 + M4
                C22_row[j] = M0_row[j] - M1_row[j] + M2_row[j] + M5_row[j];  // C22 = M1 - M2 + M3 + M6
            }
        }

        free(sub_matrices_flat);
        free(M_flat);
    }
}


int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: %s <N> <n_threads>\n", argv[0]);
        return -1;
    }
    int N = atoi(argv[1]);  // Matrix size
    int n_threads = atoi(argv[2]);

    if (N % 2 != 0) {
        printf("Matrix size must be even.\n");
        return -1;
    }

    // Allocate memory for A, B, and C
    double *A = malloc((N * N) * sizeof(double));
    double *B = malloc((N * N) * sizeof(double));
    double *C = malloc((N * N) * sizeof(double));

    // Populate A and B with random values
    populate_matrix(A, N);
    populate_matrix(B, N);

    // Print A and B matrix
    int PRINT_MATRICES = 0;
    if (PRINT_MATRICES) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%f\t", A[i * N + j]);
            }
            printf("\n");
        }
        printf("\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%f\t", B[i * N + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    double start_time = get_wall_seconds();
#pragma omp parallel num_threads(n_threads)
    {
#pragma omp single
        {
            matmul(A, B, C, N, N);
        }
    }
    double end_time = get_wall_seconds();
    printf("Matrix multiplication took %lf seconds for matrices of size %d x %d\n", end_time - start_time, N, N);

    if (PRINT_MATRICES) {
        // Print C matrix
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%f\t", C[i * N + j]);
            }
            printf("\n");
        }
    }

    free(A);
    free(B);
    free(C);
    return 0;
}