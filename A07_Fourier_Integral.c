#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cvc_numerics.h"


// Globale Konstante omega
const int w = 4;


// Fourierintegrand Cosinus für x und Parameterarray mit omega und k
double fourier_integrand_cosine(double x, void *params) {
    double w = ((double*) params)[0];
    double k = ((double*) params)[1]; 
    return cos(w*x) * cos(k*x);
}


int main(void) {
    // Integrationsparameter            
    double M[] = {10, 50, 100, 500}, left, right;           // +/- Integralgrenze M
    double delta_x = 10e-3;                                 // Weite der Integrationsschritte
    double integral_sum;                                    // Integralwert
    double parameters[2] = {w};                             // Parameter-Array für die Integrandenfunktion (Position 0: omega)
    int N_steps = 10000;                                    // Anzahl der Schritte im Intervall [-10, 10] 

    // Ausgabedatei mit Header
    FILE* fourier_file = fopen("A07_Fourier_Integral.csv", "w");
    fprintf(fourier_file, "k, Integralwert für M = 10, 30, 100, 50\n");
    
    // Itieren über die 1000 Werte von k = [-10, 10]
    double k = -10;
    for (int i = 0; i < N_steps; i++) {
        parameters[1] = k;                                  // Hinzufügen des aktuellen k zum Parameter-Arrays
        fprintf(fourier_file, "%g", k);                     // jeweiliges k in erste Spalte

        // Itieren über die 4 verschiedenen Integralgrenzen pro k
        for (int i = 0; i < 4; i++) {
            left = -M[i];                                   // untere Integralgrenze
            right = M[i];                                   // obere Integralgrenze

            // Hinzufügen der 4 Integralwerte für je ein k 
            integral_sum = cvc_integrate_simpson_2_param(left, right, delta_x, fourier_integrand_cosine, parameters);
            fprintf(fourier_file, ", %g", integral_sum);
        }
        k += 20.0 / (N_steps-1);                            // Erhöhungs des k
        fprintf(fourier_file, "\n");                        // und newline in der .csv
    }

    fclose(fourier_file);
    return 0;
}