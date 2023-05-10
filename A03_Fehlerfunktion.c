#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


#define M_PI 3.14159265358979323846


// e^(-y^2)
double e_y2(double y) {
    return 2 * exp(-pow(y, 2)) / sqrt(M_PI);
}


double cosh_midpoint(double x, double delta_x) {
    int N = fabs(x) / delta_x;
    double sum = 0, d_x = x / N;
    for (int i = 0; i < N; i++) {
        sum += cosh( (i+0.5)*d_x ) * d_x;
    }
    return sum;
}


double cosh_simpson(double x, double delta_x) {
    int N = fabs(x) / delta_x;
    double sum = 0, d_x = x / N, x_i, m_i, x_ii;
    for (int i = 0; i < N; i++) {
        x_i = d_x * i;
        m_i = d_x * (i+0.5);
        x_ii = d_x * (i+1);
        sum += ( cosh(x_i) + 4*cosh(m_i) + cosh(x_ii) ) * d_x / 6;
    }
    return sum;
}


int main(void) {
    double x, delta_x = 1 * pow(10, -4), erf_integrated_simpson, erf_integrated_midpoint;
    FILE* file_erf_simpson = fopen("data/A03_Fehlerfunktion_Simpson.csv", "w");
    fprintf(file_erf_simpson, "Integralgrenze x, Ergebnis S_erf\n");
    for (int i = 0; i < 100; i++) {
        x = -2 + i / (99.0 / 4);
        erf_integrated_simpson = erf_simpson(x, delta_x);
        erf_integrated_midpoint = erf_midpoint(x, delta_x);
        fprintf(file_erf_simpson, "%.8f, %.8f, %.8f, %.8f\n", x, erf_integrated_simpson, e_y2(x), erf_integrated_midpoint);
    }
    fclose(file_erf_simpson);

    double I = sinh(1) - sinh(0), cosh_integrated_midpoint, cosh_integrated_simpson;
    FILE* file_cosh_numint = fopen("data/A03_cosh_integrated.csv", "w");
    fprintf(file_erf_simpson, "Intervallbreite delta_x, fabs(I - S), fabs(I - M)\n");
    for (int exp = 0; exp < 10; exp++) {
        delta_x = pow(2, -exp);
        cosh_integrated_midpoint = cosh_midpoint(1, delta_x);
        cosh_integrated_simpson = cosh_simpson(1, delta_x);
        fprintf(file_cosh_numint, "%.8f, %.8f, %.8f\n", delta_x, fabs(I - cosh_integrated_simpson), fabs(I - cosh_integrated_midpoint));
    }
    fclose(file_cosh_numint);
    return 0;
}