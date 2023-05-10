#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double f2(double x) {
    return pow(x, 2) * (x - 1);
}


double diff(double x, double delta, double func(double)) {
    return ( func(x + delta) - func(x) ) / delta;
}


double centered_diff(double x, double delta, double func(double)) {
    return ( func(x + delta) - func(x - delta) ) / ( 2 * delta);
}


int main(void) {
    double f2_diff_numerical, delta, f2_centered_diff_numerical;

    FILE* file_diff_f2 = fopen("data/A05_Differenzierung.csv", "w");
    fprintf(file_diff_f2, "Delta, Abweichung\n");
    for (int exp = 0; exp < 17; exp++) {
        delta = pow(10, -exp);
        f2_diff_numerical = diff(1, delta, f2);
        f2_centered_diff_numerical = centered_diff(1, delta, f2);
        fprintf(file_diff_f2, "%.16f, %.10f, %.12f\n", delta, fabs(f2_diff_numerical - 1), fabs(f2_centered_diff_numerical - 1));
        // printf("%.16f\n", delta);
    }
    fclose(file_diff_f2);
    return 0;
}