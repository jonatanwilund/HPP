/**
 * @file strassen.h
 * @brief Utilities and public API for Strassen matrix multiplication (skeleton).
 *
 * Keep implementation details in strassen.c. Document public contracts here.
 */

#ifndef HPP_STRASSEN_H
#define HPP_STRASSEN_H

typedef struct matrix {
    double *data;  // Contiguous memory in row-major order
    int N;  // Matrix dimension (NxN)
    int tda;  // Leading dimension (stride) for each row
} matrix;

// Matrix element accessor macro
#define matrix(m, i, j) m.data[(i) * m.tda + j]

/**
 * Get wall-clock time in seconds (monotonic). Used for benchmarking.
 *
 * @return current time in seconds as a double.
 */
double get_wall_seconds(void);

/**
 * Populate a flat N-by-N matrix with test data (row-major order).
 *
 * @param A - matrix to populate
 */
void populate_matrix(matrix A);

/**
 * Perform matrix addition.
 *
 * @param A         - first matrix in addition, passed as a struct
 * @param B         - second matrix in addition, passed as a struct
 * @param C         - output matrix struct to store the result
 */
void matrix_add(const matrix *restrict A, const matrix *restrict B, matrix *restrict C);

/**
 * Perform matrix subtraction.
 *
 * @param A  - first matrix in subtraction, passed as a struct
 * @param B  - second matrix in subtraction, passed as a struct
 * @param C  - output matrix struct to store the result
 */
void matrix_subtract(const matrix *restrict A, const matrix *restrict B, matrix *restrict C);

/**
 * Perform naive matrix multiplication.
 *
 * @param A  - first matrix struct to multiply
 * @param B  - second matrix struct to multiply
 * @param C  - output matrix struct to store the result
 */
void naive_multiply(const matrix *restrict A, const matrix *restrict B, matrix *restrict C);

/**
 * Perform Strassen matrix multiplication (recursively).
 *
*
 * @param A       - first matrix struct to multiply
 * @param B       - second matrix struct to multiply
 * @param C       - output matrix struct to store the result
 */
void strassen_multiply(const matrix *restrict A, const matrix *restrict B, matrix *restrict C);

#endif // HPP_STRASSEN_H
