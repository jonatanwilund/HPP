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
    unsigned int N;  // Matrix dimension (NxN)
    unsigned int tda;  // Leading dimension (stride) for each row
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
 * @param C         - output matrix descriptor
 * @param A         - first matrix descriptor
 * @param B         - second matrix descriptor
 */
void matrix_add(matrix C, matrix A, matrix B);

/**
 * Perform matrix subtraction.
 *
 * @param C  - output matrix descriptor
 * @param A  - first matrix descriptor
 * @param B  - second matrix descriptor
 */
void matrix_subtract(matrix C, matrix A, matrix B);

/**
 * Assemble matrix from sub-matrices.
 *
 * @param C     - output matrix struct to store the result
 * @param C11   - top-left sub-matrix
 * @param C12   - top-right sub-matrix
 * @param C21   - bottom-left sub-matrix
 * @param C22   - bottom-right sub-matrix
 */
void assemble_matrix(
    matrix C,
    matrix C11,
    matrix C12,
    matrix C21,
    matrix C22);

/**
 * Perform naive matrix multiplication.
 *
 * @param C  - output matrix descriptor
 * @param A  - first matrix descriptor
 * @param B  - second matrix descriptor
 */
void naive_multiply(matrix C, matrix A, matrix B);

/**
 * Perform Strassen matrix multiplication (recursively).
 *
*
 * @param C       - output matrix descriptor
 * @param A       - first matrix descriptor
 * @param B       - second matrix descriptor
 */
void strassen_multiply(matrix C, matrix A, matrix B);

#endif // HPP_STRASSEN_H
