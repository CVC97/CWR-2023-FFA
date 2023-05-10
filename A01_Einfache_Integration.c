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


// numerische Integration
double integrate(double x, int N) {
    double sum = 0, interval = x / N;
    for (int i = 0; i < N; i++) {
        sum += p1(interval * i) * interval;
    }
    return sum;
}


// Aufruf der numerischen Integration
int main(void) {
    int N = 1000;          // Anzahl der Teilintervalle
    double x = 8;      // Intervallgrenze
    double analytical_solution = p1_integrated(x);
    double numerical_solution = integrate(x, N); 
    printf("Analytische Lösung P1(x = %.2f):\t %.2f\n", x, analytical_solution);
    printf("Numerische Lösung (N = %d, x = %.2f):\t %.2f\n", N, x, numerical_solution);
    return 0;
}