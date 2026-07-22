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
 * Perform naive matrix multiplication.
 *
 * A, B    - input matrices in row-major order (size n*n)
 * C       - output matrix buffer (size n*n), owned by caller
 * N       - matrix dimension (must be a power of two)
 *
 * @return 0 on success, non-zero on error.
 */
int naive_multiply(const double *A, const double *B, double *C, unsigned int N);

/**
 * Perform Strassen matrix multiplication.
 *
 * A, B    - input matrices in row-major order (size n*n)
 * C       - output matrix buffer (size n*n), owned by caller
 * n       - matrix dimension (must be a power of two)
 *
 * @return 0 on success, non-zero on error.
 */
int strassen_multiply(const double *A, const double *B, double *C, unsigned int N);

#endif // HPP_STRASSEN_H
