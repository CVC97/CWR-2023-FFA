#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cvc_numerics.h"


int main(void) {
    int N = 100000;
    long int length = 1000000;

    // malloc
    time_t t_malloc_0 = time(0);
    for (int i = 0; i < N; i++) {
        double *malloc_array = (double*) malloc(length * sizeof(double));
        if (malloc_array == NULL) {
            printf("ERROR! RAM FULL.\n");
            return 1;
        }
        for (int j = 0; j < length; j++) {
            malloc_array[j] = 0;
        }
        free(malloc_array);
    }
    time_t t_malloc_1 = time(0);
    printf("malloc %d x %ld: %li\n", N, length, t_malloc_1 - t_malloc_0);

    // calloc
    time_t t_calloc_0 = time(0);
    for (int i = 0; i < N; i++) {
        double *calloc_array = (double*) calloc(length, sizeof(double));
        if (calloc_array == NULL) {
            printf("ERROR! RAM FULL.\n");
            return 1;
        }
        free(calloc_array);
    }
    time_t t_calloc_1 = time(0);
    printf("calloc %d x %ld: %li\n", N, length, t_calloc_1 - t_calloc_0);

    return 0;
}