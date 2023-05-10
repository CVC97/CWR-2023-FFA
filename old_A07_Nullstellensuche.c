#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double f(double x, double t) {
    return x - tanh(x/t);
}


// numerische Differenzierung: zentrale Differenz
double centered_diff(double x, double t, double delta, double func(double, double)) {
    return ( func(x + delta, t) - func(x - delta, t) ) / ( 2 * delta);
}


// Nullstellensuche: Newton-Raphson
double find_root_newt_raph(double func(double, double), double t, double x0, double delta, double rel_tol, int max_iter) {
    int i = 0;
    double x_old;

    FILE* file_root_newt_raph = fopen("data/old_A07_Nullstelle_Newton_Raphson.csv", "w");
    fprintf(file_root_newt_raph, "Itierungsschritt, abs(x_(n+1) - x_n)\n");
    while (i++ < max_iter) {
        printf("i: %d, max_iter: %d\n", i, max_iter);
        x_old = x0;
        x0 -= func(x0, t) / centered_diff(x0, t, delta, func);
        fprintf(file_root_newt_raph, "%d, %.6f\n", i, fabs(x0 - x_old));
        printf("%d, %.6f\n", i, fabs(x0 - x_old));
        if (fabs(x_old - x0) / fabs(x0) < rel_tol) {
            break;
        }
    }
    fclose(file_root_newt_raph);
    return x0;
}


// Nullstellensuche: Newton-Raphson
double find_root_newt_raph_t(double func(double, double), double t, double x0, double delta, double rel_tol, int max_iter) {
    int i = 0;
    double x_old;
    while (i++ < max_iter) {
        x_old = x0;
        x0 -= func(x0, t) / centered_diff(x0, t, delta, func);
        if (fabs(x_old - x0) / fabs(x0) < rel_tol) {
            break;
        }
    }
    return x0;
}


// Nullstellensuche: Bisektion
double find_root_bisection(double func(double, double), double a, double b, double epsilon, int max_iter) {
    double t = 0.5;
    double x_mid, f_x_mid, f_a, f_b, x_mid_old;
    int i = 0;

    FILE* file_root_bisection = fopen("data/old_A07_Nullstelle_Bisektion.csv", "w");
    fprintf(file_root_bisection, "Itierungsschritt, abs(x_(n+1) - x_n)\n");
    while (i++ < max_iter) {
        x_mid_old = x_mid;
        x_mid = (a + b) / 2;
        f_x_mid = func(x_mid, t);
        f_a = func(a, t);
        f_b = func(b, t);
        fprintf(file_root_bisection, "%d, %.6f\n", i, fabs(x_mid - x_mid_old));
        if (f_x_mid * f_b < 0) {
            a = x_mid;
        } else if (f_x_mid * f_a < 0) {
            b = x_mid;
        }
        if (f_a < epsilon && f_b < epsilon) {
            break;
        }
    }
    fclose(file_root_bisection);
    return x_mid;
}


int main(void) {
    double root_newt_raph, root_bisection, root_t_newt_raph;

    root_newt_raph = find_root_newt_raph(f, 0.5, 1, 0.1, 0.00001, 20);
    root_bisection = find_root_bisection(f, 0.1, 4, 0.00001, 20);

    // Ausgabe der Nullstelle auf 6 Ziffern nach dem Komma
    printf("Nullstelle (Newton-Raphson): %.6f\n", root_newt_raph);
    printf("Nullstelle (Bisektion): %.6f\n", root_bisection);
    
    // Abhängigkeit der Nullstelle von t (Newton-Raphson, t in [0, 1.2])
    FILE* file_root_t_newt_raph_t = fopen("data/old_A07_Nullstelle_t_Newton_Raphson.csv", "w");
    fprintf(file_root_t_newt_raph_t, "t, Nullstelle\n");
    for (double t = 0; t <= 1.2; t += 0.01) {
        root_t_newt_raph = find_root_newt_raph_t(f, t, 1, 0.1, 0.00001, 100);
        fprintf(file_root_t_newt_raph_t, "%f, %.6f\n", t, root_t_newt_raph);
    }
    fclose(file_root_t_newt_raph_t);
    return 0;
}