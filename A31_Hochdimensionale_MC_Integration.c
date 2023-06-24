#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "../cvc_numerics.h"
#include "../cvc_rng.h"


// Unity Sphere density in D dimensions
double density_func(int D, double x[]) {
    double sum_x = 0;
    for (int i_dim = 0; i_dim < D; i_dim++) {
        sum_x += cvc_npow(x[i_dim], 2);
    }
    if (sum_x > 1) {
        return 0;
    } else {
        return 1;
    }
}


// Numerisches Volumen der D-dimensionalen Kugel mit Radius r
double volume_unity_sphere_analytical(int D, double r) {
    double Gamma;
    if (D % 2 == 0) {
        Gamma = cvc_factorial(D/2);
    } else {
        double n = D/2.0 + 0.5;
        Gamma = sqrt(cvc_PI) * cvc_factorial(2*n) / cvc_factorial(n) / cvc_npow(4, n);
    }
    return cvc_npow(r, D) * pow(cvc_PI, D/2.0) / Gamma;
}


int main(void) {

    // +++++ (Aufgabe 4: Numerische MC-Integration und analytische Lösung des Volumens der Einheitssphäre in 1-7 Dimensionen) +++++
    {
    
    int N = 10e5;                                                                       // Anzahl der Stützstellen für die MC-Integration
    FILE* volume_dim_file = fopen("A31_volume_dimension.csv", "w");                     // Datenfile
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);                                      // Wahl des Zufallsgenerators: gsl_rng_mt19937

    // Iteration über Dimensionen 1-7
    for (int dimension = 1; dimension < 8; dimension++) {                               // Anzahl Dimensionen
        printf("calculating volume %d/7...\n", dimension);
        fprintf(volume_dim_file, "%d", dimension);                             

        // Maße des Hyperquaders [-1, 1] für D-Dimensionen
        double R[2 * dimension];
        for (int i_dim = 0; i_dim < dimension; i_dim++) {
            R[2*i_dim] = -1;
            R[2*i_dim + 1] = 1;
        }

        // 100-malige Wiederholung der MC-Integration für jede Dimension
        for (int i = 0; i < 10; i++) {
            double numerical_volume = cvc_mc_integrate(rng, dimension, density_func, R, N);
            fprintf(volume_dim_file, ", %g", numerical_volume);
        }
        double analytical_volume = volume_unity_sphere_analytical(dimension, 1);
        fprintf(volume_dim_file, ", %g\n", analytical_volume);  
    }
    fclose(volume_dim_file);
    free(rng);


    }
    // +++++ (Aufgabe 5: Differenz der numerische MC-Integration zur analytische Lösung des Volumens der Einheitssphäre in 5 Dimensionen für verschiedene Anzahlen von Stützstellen) +++++
    {

        int dimension = 5;                                                              // Anzahl der Dimensionen
        double N_exp_min = 2, N_exp_max = 8;                                            // Laufvariablen der Stützstellen für die MC-Integration
        double analytical_volume = volume_unity_sphere_analytical(dimension, 1);        // Analytisches Ergebnis

        FILE* residue_n_file = fopen("A31_residue_n.csv", "w");                         // Datenfile
        fprintf(residue_n_file, "N Stützstellen, Abweichung zum analytischen Ergebnis\n");
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);                                  // Wahl des Zufallsgenerators: gsl_rng_mt19937

        // Maße des Hyperquaders [-1, 1] für 5 Dimensionen
        double R[2 * dimension];
        for (int i_dim = 0; i_dim < dimension; i_dim++) {
            R[2*i_dim] = -1;
            R[2*i_dim + 1] = 1;
        }

        // Iteration über 100 normalverteilte Anzahlen an Stützstelen von 10e2 bis 10e8
        for (int i_exp = 0; i_exp < 100; i_exp++) {
            printf("calculating residue %d/100...\n", i_exp + 1);
            double N_exp = N_exp_min + (N_exp_max - N_exp_min) * i_exp / 99.0;
            int N = pow(10, N_exp);
               
            double numerical_volume = cvc_mc_integrate(rng, dimension, density_func, R, N);
            fprintf(residue_n_file, "%d, %g\n", N, fabs(numerical_volume - analytical_volume));
        }
        fclose(residue_n_file);
        free(rng);


    }
    // +++++ (Aufgabe 6: Standardabweichung des Volumens der Einheitssphäre in 5 Dimensionen für verschiedene Anzahlen von Stützstellen und M = 100 Datensätze) +++++
    {

        int dimension = 5;                                                              // Anzahl der Dimensionen
        double N_exp_min = 2, N_exp_max = 8;                                            // Laufvariablen der Stützstellen für die MC-Integration
        int M = 100;                                                                    // Anzahl der Datenpunkte pro Stützstellenzahl N
        double volume_array[M];                                                         // Array zur temporären Speicherung der Datenpunkte

        FILE* sigma_n_file = fopen("A31_sigma_n.csv", "w");                             // Datenfile
        fprintf(sigma_n_file, "N Stützstellen, Standardabweichung\n");
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);                                  // Wahl des Zufallsgenerators: gsl_rng_mt19937

        // Maße des Hyperquaders [-1, 1] für 5 Dimensionen
        double R[2 * dimension];
        for (int i_dim = 0; i_dim < dimension; i_dim++) {
            R[2*i_dim] = -1;
            R[2*i_dim + 1] = 1;
        }

        // Iteration über 100 normalverteilte Anzahlen an Stützstelen von 10 bis 10e6
        for (int i_exp = 0; i_exp < 100; i_exp++) {
            printf("calculating sigma %d/100...\n", i_exp + 1);
            double N_exp = N_exp_min + (N_exp_max - N_exp_min) * i_exp / 99.0;
            int N = pow(10, N_exp);
            double volume_sum = 0;

            // Erzeugung M = 100 Datenpunkte pro N
            for (int i_M = 0; i_M < M; i_M++) {
                double numerical_volume = cvc_mc_integrate(rng, dimension, density_func, R, N / M);
                volume_sum += numerical_volume;
                volume_array[i_M] = numerical_volume;
            }

            // Berechnung der Standardabweichung mit den M = 100 Arrayeinträgen
            double volume_mean = volume_sum / M;                                        // Mittelwert der 100 Datenpunkte
            double sigma_square_sum = 0;                                                // Berechnung der Standardabweichung
            for (int i_M = 0; i_M < M; i_M++) {                                         // Erzeugung M = 100 Datensätze
                sigma_square_sum += cvc_npow(volume_array[i_M] - volume_mean, 2);
            }
            fprintf(sigma_n_file, "%d, %g\n", N, sqrt(sigma_square_sum / M));
        }
        fclose(sigma_n_file);
        free(rng);


    }
    return 0;
}