#include <sys/time.h>
#include "strassen.h"


double get_wall_seconds() {
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