#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Funktion p1(x)
double p1(double x) {
    return pow(x, 3) - x / 2;
}


// analytische Lösung des Integrals von p_1(x)
double p1_integrated(double x) {
    return pow(x, 4) / 4 - pow(x, 2) / 4;
}


// numerische Integration (Untersumme)
double integrate(double left, double right, int N, double integrand(double)) {
    double sum = 0, interval = (right - left) / N;
    for (int i = 0; i < N; i++) {
        sum += integrand(left + i*interval) * interval;
    }
    return sum;
}


// numerische Integration (Trapez)
double integrate_pez(double left, double right, int N, double integrand(double)) {
    double sum = 0, interval = (right - left) / N;
    for (int i = 0; i < N; i++) {
        sum += (integrand(left + i*interval) + integrand(left + (i+1)*interval)) * interval / 2;
    }
    return sum;
}


// Aufruf der numerischen Integration und Erzeugen der Datenfiles zu Auswertung
int main(void) {
    double left = 0, right = 2;      
    double analytical_solution = p1_integrated(right) - p1_integrated(left), numerical_solution_left, numerical_solution_pez;

    printf("Intervall: [%.4f, %.4f]\n", left, right);
    printf("Analytische Lösung:\t\t %.4f\n", analytical_solution);

    FILE* file_numerical_left = fopen("data/A02_Untersumme.csv", "w");
    FILE* file_numerical_pez = fopen("data/A02_Trapez.csv", "w");

    fprintf(file_numerical_left, "Teilintervalle N, Integral [%.4f, %.4f]\n", left, right);
    fprintf(file_numerical_pez, "Teilintervalle N, Integral [%.4f, %.4f]\n", left, right);

    int N = 1;
    for (int i = 0; i < 7; i++) {
        N *= 10;
        numerical_solution_left = integrate(left, right, N, p1); 
        numerical_solution_pez = integrate_pez(left, right, N, p1); 

        printf("Numerische Lösung (N = 10^%d):\t %.4f (Untersumme)\t %.4f (Trapez)\n", i+1, numerical_solution_left, numerical_solution_pez);

        fprintf(file_numerical_left, "%d, %.8f\n", N, numerical_solution_left - analytical_solution);
        fprintf(file_numerical_pez, "%d, %.8f\n", N, numerical_solution_pez - analytical_solution);
    }
    fclose(file_numerical_left);
    fclose(file_numerical_pez);
    return 0;
}