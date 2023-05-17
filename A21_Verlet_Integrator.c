#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"

typedef int ode_func(double, const double[], double[], void*);


/*------------------------  PHYSIKALISCHE KONSTANTEN  -------------------------*/

const int N = 30;                       // Anzahl der Pendel
const double k = 100;                   // Federhaerte [N/m]
const double base_length = 0;           // Basislaenge der Federn [m]
const double mass = 1;                  // Masse der Pendel [kg]

/*-------------------------  SIMULATIONS-PARAMETER  ---------------------------*/

const double T_max = 20.0;
const double delta_t = 1e-3; // Zeitliche Schrittweite


// Lösen der ODE
int pendumlumsODE(double t, const double y[], double f[], void *params) {
    // Geschwindigkeiten aus dem Zustandsarray übertragen
    for (int i = 0; i < N; i++) {                   
        f[i] = y[N+i];
    }

    // Beschleunigungen aus Positionen berechnen
    f[N] = -k*y[0] / mass;                                      // Position 0 generell
    if (N != 1) {
        f[N] += k*(y[1] - y[0] - 1);                            // Position 1 auf Position 0 wenn mehr als 1 Pendel
        f[2*N-1] += -k*(y[N-1] - y[N-2] - 1) / mass;            // Position N - 1 (wenn mehr als 1 Pendel)
    }
    for (int i = 1; i < N - 1; i++) {                           // Positionen 1 bis N - 2
        f[N+i] += -k*(y[i] - y[i-1] - 1) / mass + k*(y[i+1] - y[i] - 1) / mass; 
    }
    return 0;
}


// Berechnung der Energie
double pendulums_energy(const double y[]) {
    double E_pot = 0;
    double E_kin = 0;

    // Kinetische Energie berechnen
    for (int i = 0; i < N; i++) {
        E_kin += mass / 2 * npow(y[N+i], 2);
    }

    // Potentielle Energie
    E_pot += k / 2 * npow(y[0], 2);  
    if (N != 1) {
        E_pot += k / 2 * npow(y[1] - y[0] - 1, 2);
        E_pot += k / 2 * npow(y[N-1] - y[N-2] - 1, 2);                                     
    }
    for (int i = 1; i < N - 1; i++) {                          
        E_pot += k / 2 * npow(y[i] - y[i-1] - 1, 2) + k / 2 * npow(y[i+1] - y[i] - 1, 2); 
    }
    return E_kin + E_pot;
}


int main(void) {
    int dimension = 2*N;
    double t = 0;

    // Erstellen des Zustandsarrays mit den Startwerten
    double y[dimension], energy;
    for (int i = 0; i < N; i++) {
        y[i] = i;
        y[N+i] = 0;
    }
    y[N] = 20;

    // Ausgabe-Dateien
    FILE* pos_file = fopen("data/A21_Verlet_Integrator.csv", "w");
    FILE* energy_file = fopen("data/A21_Verlet_Integrator_Energie.csv", "w");
    fprintf(pos_file, "Zeit t, Pendel");
    fprintf(energy_file, "Zeit t, Pendel");
    for (int i = 1; i <= N; i++) {
        fprintf(pos_file, ", P%d", i);
        fprintf(energy_file, ", P%d", i);
    }
    fprintf(pos_file, "\n");
    fprintf(energy_file, "\n");


    // Durchlaufen der Zeitschritte
    while (t < T_max) {
        verlet_step(t, delta_t, y, pendumlumsODE, dimension, NULL);     // Aufruf der Integrationsfunktion
        energy = pendulums_energy( (const double*) y);                  // Berechnung der GEsamtenergie des Systems

        // Beschreiben des Datenfiles
        fprintf(pos_file, "%g", t);
        fprintf(energy_file, "%g, %g\n", t, energy);
        for (int i = 0; i < N; i++) {
            fprintf(pos_file, ", %g", y[i]);
        }
        fprintf(pos_file, "\n"); 

        t += delta_t;
    }

    fclose(pos_file);
    fclose(energy_file);
    return 0;
}