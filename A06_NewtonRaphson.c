#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


const double mu_earth = 3.986e14;
const double mu_sun = 1.327e20;
const double omega = 1.991e-7;
const double radius_orbit_earth = 1.496e11;


double acceleration(double r) {
    double a_earth = -mu_earth / pow(r - radius_orbit_earth, 2);
    double a_sun = -mu_sun / pow(r, 2);
    double a_centrifugal = pow(omega, 2) * r;
    return a_earth + a_sun + a_centrifugal;
}


int main(void) {
    double lagrange, x0 = 1.5e11, delta = 1000000, rel_tol = 0.00000001, max_iter = 1000;
    lagrange = find_root_newton_raphson(acceleration, x0, delta, rel_tol, max_iter) - radius_orbit_earth;
    printf("Lagrangepunkt L2: %.2f Mio. km \n", lagrange / 1000000000);
    return 0;
}