/*---------------------------------------------*\
|  CWR SoSe 2022                                     |
|  Blatt 5 - Aufgabe 15                         |
|  Contact: salvatore.manmana@uni-goettingen.de   |
|  Original Author: mail@saschalambert.de       |
\*---------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "cvc_numerics.h"

/* TODO */
/*------------------------  PHYSIKALISCHE KONSTANTEN  -------------------------*/

const int N = 30;                       // Anzahl der Pendel
const double k = 100;                   // Federhaerte [N/m]
const double base_length = 0;           // Basislaenge der Federn [m]
const double mass = 1;                  // Masse der Pendel [kg]

/*-------------------------  SIMULATIONS-PARAMETER  ---------------------------*/

const double T_max = 20.0;
const double delta_t = 1e-2; // Zeitliche Schrittweite
const double smoothing_factor = 10e-10;

static const char pos_file_name[] = "data/A12_position.csv";
static const char energy_file_name[] = "data/A12_energy.csv";

/*-------------------------  PHYSIKALISCHES SYSTEM  ---------------------------*/

/* TODO */
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


/* TODO */
/* Bearbeiten Sie hier Aufgabe 15.9 */
double pendulums_energy(const double y[]) {
    double E_pot = 0;
    double E_kin = 0;

    // Kinetische Energie berechnen
    for (int i = 0; i < N; i++) {
        E_kin += mass / 2 * cvc_npow(y[N+i], 2);
    }

    // Potentielle Energie
    E_pot += k / 2 * cvc_npow(y[0], 2);  
    if (N != 1) {
        E_pot += k / 2 * cvc_npow(y[1] - y[0] - 1, 2);
        E_pot += k / 2 * cvc_npow(y[N-1] - y[N-2] - 1, 2);                                     
    }
    for (int i = 1; i < N - 1; i++) {                          
        E_pot += k / 2 * cvc_npow(y[i] - y[i-1] - 1, 2) + k / 2 * cvc_npow(y[i+1] - y[i] - 1, 2); 
    }
    return E_kin + E_pot;
}

/*-----------------------------------------------------------------------------*/

int main(void)
{
    /* ------- GSL Verwaltung -------------------------------------------------*/

    /* TODO */
    /* Berechnen Sie hier die korrekte Dimensionalitaet des DGL-Systems */
    int gsl_dimension = 2*N;

    /* Initialisieren des Systems fuer GSL
       Wir uebergeben Ihre erstellte Funktion und die Dimensionaliaet */
    gsl_odeiv2_system SYSTEM = {pendumlumsODE, NULL, gsl_dimension, NULL};

    /* Auswahl des Integrators
       Wir waehlen hier einen Runge-Kutta 4. Ordnung (rk4). */
    gsl_odeiv2_step *STEPPER = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, gsl_dimension);

    /* ------- Anfangswerte ---------------------------------------------------*/

    /* TODO */
    /* Erstellen Sie zwei Arrays mit der notwendigen Groesse. Sie muessen
       alle Positionen und Geschwindigkeiten fassen koennen. */
    double y[gsl_dimension];
    double yerr[gsl_dimension]; /* GSL benoetigt yerr */

    /* TODO */
    /* Fuellen Sie das Array y[] mit den Startwerten */
    // double y[] = {
    //     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    //     20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
    // };
    for (int i = 0; i < N; i++) {
        y[i] = i;
        y[N+i] = 0;
    }
    y[N] = 20;

    /* ------- Simulations-Schleife -------------------------------------------*/

    /* Ausgabe-Dateien oeffnen */
    FILE *pos_file = fopen(pos_file_name, "w");
    FILE *energy_file = fopen(energy_file_name, "w");


    // Header der Ausgabe-Datei
    fprintf(pos_file, "Zeit t");
    fprintf(energy_file, "Teit t, Energie E\n");
    for (int i = 1; i <= N; i++) {
        fprintf(pos_file, ", P%d", i);
    }
    fprintf(pos_file, "\n"); 

    /* Simulationszeit */
    double t = 0;

    /* #####   HAUPTSCHLEIFE   ##### */
    while (t < T_max)
    {
        /* step_apply befoertdert y[] zum naechsten Zeitpunkt */
        gsl_odeiv2_step_apply(STEPPER, t, delta_t, y, yerr, NULL, NULL, &SYSTEM);

        // Berechnung der Energie
        double energy = pendulums_energy( (const double*) y);
        

        /* TODO */
        /* Geben Sie hier die Daten in die Dateien aus */
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
    return EXIT_SUCCESS;
}
