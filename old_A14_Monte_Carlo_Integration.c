#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "cvc_numerics.h"


// Unity Sphere density
double rho_Sn(double y[], int dim) {
    double sum_y = 0;
    for (int i = 0; i < dim; i++) {
        sum_y += pow(y[i], 2);
    }
    if (sum_y > 1) {
        return 0;
    }
    return 1;
}


// Radial Symmetry density
double rho_Rad(double y[], int dim) {
    double sum_y = 0;
    for (int i = 0; i < dim; i++) {
        sum_y += pow(y[i], 2);
    }
    if (sum_y > 1 || sum_y <= 0.9) {
        return 0;
    }
    return 1;
}


int main(void) {
    double a[] = {-8, -1, -1, -4, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    double b[] = {2, 4, 2, 5, 10, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, -1, -1, -1, -1, -1};
    double volume_us, volume_rs;

    // Unity Sphere (||y|| <= 1)
    FILE* file_unity_sphere_volume = fopen("data/A14_Volumen_Einheitskugel.csv","w");
    fprintf(file_unity_sphere_volume, "Dimensionen, Volumen\n");
    for (int dim = 1; dim <= 15; dim++) {
        volume_us = mc_integrate(rho_Sn, a, b, dim, 1e6);
        fprintf(file_unity_sphere_volume, "%d, %f\n", dim, volume_us);
        // printf("Volumen Einheitsspäre (%d Dimensionen): %f\n", dim, volume_us);
    }
    fclose(file_unity_sphere_volume);

    // Radial Symmetry (0.9 < ||y|| <= 1)
    FILE* file_radial_symmetry_volume = fopen("data/A14_Volumen_Radialsymmetrie.csv","w");
    fprintf(file_radial_symmetry_volume, "Dimensionen, Volumen\n");
    for (int dim = 1; dim <= 18; dim++) {
        volume_rs = mc_integrate(rho_Rad, a, b, dim, 1e8);
        fprintf(file_radial_symmetry_volume, "%d, %f\n", dim, volume_rs);
        // printf("Volumen Radialsymmetrie (%d Dimensionen): %f\n", dim, volume_rs);
    }
    fclose(file_radial_symmetry_volume);
    return 0;
}