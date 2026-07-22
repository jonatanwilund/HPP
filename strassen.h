/**
 * @file strassen.h
 * @brief Utilities and public API for Strassen matrix multiplication (skeleton).
 *
 * Keep implementation details in strassen.c. Document public contracts here.
 */

#ifndef HPP_STRASSEN_H
#define HPP_STRASSEN_H

/**
 * Get wall-clock time in seconds (monotonic). Used for benchmarking.
 *
 * @return current time in seconds as a double.
 */
double get_wall_seconds(void);

/**
 * Populate a flat N-by-N matrix with test data (row-major order).
 *
 * @param A buffer of at least N*N doubles (owned by caller)
 * @param N matrix dimension
 */
void populate_matrix(double *A, int N);

/**
 * Perform matrix addition.
 *
 * @param A         - input matrix in row-major order (size NxN)
 * @param B         - input matrix in row-major order (size NxN)
 * @param C         - output matrix buffer (size NxN), owned by caller
 * @param N         - matrix dimension
 * @param stride    - memory stride for each row (equal to the original dimension of the multiplied matrices)
 */
void matrix_add(
    const double *restrict A,
    const double *restrict B,
    double *restrict C,
    unsigned int N,
    unsigned int stride);

/**
 * Perform matrix subtraction.
 *
 * @param A  - input matrix in row-major order (size NxN)
 * @param B  - input matrix in row-major order (size NxN)
 * @param C  - output matrix buffer (size NxN), owned by caller
 * @param N  - matrix dimension
 * @param stride    - memory stride for each row (equal to the original dimension of the multiplied matrices)
 */
void matrix_subtract(
    const double *restrict A,
    const double *restrict B,
    double *restrict C,
    unsigned int N,
    unsigned int stride);

/**
 * Perform naive matrix multiplication.
*
 * @param A  - input matrix in row-major order (size NxN)
 * @param B  - input matrix in row-major order (size NxN)
 * @param C  - output matrix buffer (size NxN), owned by caller
 * @param N  - matrix dimension
 */
void naive_multiply(const double *restrict A, const double *restrict B, double *restrict C, unsigned int N);

/**
 * Perform Strassen matrix multiplication.
 *
*
 * @param A       - input matrix in row-major order (size NxN)
 * @param B       - input matrix in row-major order (size NxN)
 * @param C       - output matrix buffer (size n*n), owned by caller
 * @param N       - matrix dimension
 * @param stride  - memory stride for each row (equal to the original dimension of the multiplied matrices)
 */
void strassen_multiply(
    const double *restrict A,
    const double *restrict B,
    double *restrict C,
    unsigned int N,
    unsigned int stride);

#endif // HPP_STRASSEN_H
