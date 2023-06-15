#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gsl/gsl_rng.h>


// Struktur des inneren Zustandes eines Zufallsgenerators
struct randu_random_t {
    uint32_t state;
};


// Zufallsgenerator des Typs RANDNU
uint32_t randu_random_r(struct randu_random_t *rng) {
    a = (1u << 16) + 3;
    b = 0;
    m =  (1u << 31);
}


int main(void) {

    return 0;
}