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
 * @param depth   - current recursion depth
 */
void strassen_multiply_recursive(matrix C, matrix A, matrix B, unsigned int depth);

#endif // HPP_STRASSEN_H
