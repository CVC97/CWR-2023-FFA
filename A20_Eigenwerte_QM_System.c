#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cvc_numerics.h"


// Quantenzahlen / physikalische Konstanten
const int m = 1;
const int k = 1;
const int h_cross = 1;


// rechte Seite der zeitunabhängigen Schrödinger-Gleichung
int schroedinger_harmonic_ODE(double x, const double y[], double f[], void *params) {
    double E = *(double*) params;                                                       // Übergabe des Energieeigenwertes über den void-Pointer params
    double V = 1.0 / 2 * k * cvc_npow(x, 2);                                            // Berechnung des Potentials für gegebenes x

    f[0] = y[1];                                                                        // Übertragen der ersten Ableitung Psi-Prime in Ableitungsarray
    f[1] = -2 * m / cvc_npow(h_cross, 2) * (E - V) * y[0];                              // Berechnung der zweiten Ableitung
    return 0;
}


// Rückgabe des Produktes aus Psi(x1) * Psi-Prime(x1) für gegebenen Energiewert E
double test_energy(double E) {
    double psi[2] = {1, 0};                                                             // Initialisierung der Startbedingungen: Psi(x0) = 1, Psi_Prime(x0) = 0
    double params[1] = {E};                                                             // Initialisieren des Parameter-Arrays mit dem Energiewert

    // Numerische Integration von x0 = -1 bis x1 = 0 mit RK4
    double x0 = -10;                                                                    // Startwert der Integration
    double x1 = 0;                                                                      // Endpunkt der Integration
    double delta_x = 1e-3;                                                              // Schrittweite der Integration
    while (x0 + delta_x < x1) {
        cvc_rk4_step(x0, delta_x, psi, schroedinger_harmonic_ODE, 2, params);
        x0 += delta_x;
    }
    cvc_rk4_step(x0, x1 - x0, psi, schroedinger_harmonic_ODE, 2, params);
    return psi[0] * psi[1];                                                             // Rückgabe: Psi(x1) * Psi-Prime(x1) = 0 wenn E Energieeigenwert von Psi 
} 


int main(void) {

// +++++ (Aufgabe 3: Itierung über Energiewerte im Intervall [0.1, 10]) +++++

    double E_start = 0.1;   
    double E_end = 10;

    FILE* test_energy_file = fopen("A20_test_energy.csv", "w");
    fprintf(test_energy_file, "Energiewert, Psi(x1) * Psi_Prime(x1)\n");
    
    int N = 1000;                                                                       // Anzahl der Energiewerte im Intervall
    for (int i = 0; i < N; i++) {
        double E = E_start + (E_end - E_start) * i / (N - 1);
        fprintf(test_energy_file, "%g, %g\n", E, test_energy(E));
    }
    fclose(test_energy_file);


// +++++ (Aufgabe 4: Bestimmung der Nullpunktsenergien für n = 0 und n = 1 mit Newton-Raphson) +++++

    double E_newt_raph_0 = 0.4;                                                         // Startenergie zur Nullstellensuche des ersten Energieigenwertes E_0
    double E_newt_raph_1 = 1.4;                                                         // Startenergie zur Nullstellensuche des zweiten Energieigenwertes E_1
    
    double delta_E = 1e-3;                                                              // delta_E der numerischen Differenzierung in der Nullsellensuche
    double rel_tol = 1e-3;                                                              // Relative Toleranz zum Abbruch
    double max_iter = 1e4;                                                              // Maximale Anzahl der Schleifendurchläufe der Nullstellensuche

    double E_0 = cvc_find_root_newton_raphson(test_energy, E_newt_raph_0, delta_E, rel_tol, max_iter);
    double E_1 = cvc_find_root_newton_raphson(test_energy, E_newt_raph_1, delta_E, rel_tol, max_iter);

    printf("E_0: %.10g\n", E_0);
    printf("E_1: %.10g\n", E_1);


// +++++ (Aufgabe 4: Numerische Berechnung der Wellenfunktionen der Nullpunktsenergien E_0 und E_1 für n = 0 und n = 1) +++++

    double psi_0[2] = {1, 0};                                                           // Startbedinungen der Integration für Energieeigenwert E_0
    double psi_1[2] = {1, 0};                                                           // Startbedinungen der Integration für Energieeigenwert E_1

    double params_0[1] = {E_0};                                                         // Übergabearray mit Energieeigenwert E_0
    double params_1[1] = {E_1};                                                         // Übergabearray mit Energieeigenwert E_1

    FILE* psi_file = fopen("A20_psi_x.csv", "w");
    fprintf(psi_file, "x, Psi_0(x), Psi_1(x)\n");

    // Numerische Integration von x0 = -1 bis x1 = 0 mit Runke-Kutta-4
    double x0 = -10;                                                                    // Startwert der Integration
    double x1 = 0;                                                                      // Endpunkt der Integration
    double delta_x = 1e-3;                                                              // Schrittweite der Integration
    fprintf(psi_file, "%g, %g, %g\n", x0, psi_0[0], psi_1[0]);
    while (x0 + delta_x < x1) {
        cvc_rk4_step(x0, delta_x, psi_0, schroedinger_harmonic_ODE, 2, params_0);
        cvc_rk4_step(x0, delta_x, psi_1, schroedinger_harmonic_ODE, 2, params_1);
        x0 += delta_x;
        fprintf(psi_file, "%g, %g, %g\n", x0, psi_0[0], psi_1[0]);
    }
    cvc_rk4_step(x0, x1 - (x0 + delta_x), psi_0, schroedinger_harmonic_ODE, 2, params_0);
    cvc_rk4_step(x0, x1 - (x0 + delta_x), psi_1, schroedinger_harmonic_ODE, 2, params_1);
    fprintf(psi_file, "%g, %g, %g\n", x0, psi_0[0], psi_1[0]);
    fclose(psi_file);
    return 0;
}