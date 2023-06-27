#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "cvc_numerics.h"
#include "cvc_rng.h"


// How much is the sickness?
int how_sick_is_the_bus(int bus[], int N) {
    int sick_sum = 0;
    for (int i_N = 0; i_N < N; i_N++) {
        sick_sum += bus[i_N];
    }
    return sick_sum;
}


int main(void) {

    // Parameter des Systems
    int N = 40;                                                 // Anzahl der Sitzplätze im Bus N
    int bus_array[N];                                           // Array mit Zuständen der Passagiere
    double sick_0 = 0.2;                                        // Wahrscheinlichkeit sick_0 eines Passagieres krank zu sein bei t = 0
    double p = 0.35;                                            // Wahrscheinlichkeit einer Spontangenesung p
    int T_max = 1000;                                           // Anzahl der Zeitschritte T_max

    // Wahl des Zufallsgenerators: gsl_rng_mt19937
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL)); 

    // Datenfile und erste Zeile mit Anzahl der Zeitschritten
    FILE* bus_file = fopen("data/A32_Bus_Schnupfen.csv", "w");
    for (int i_T = 0; i_T < T_max; i_T++) {
        fprintf(bus_file, "%d, ", i_T);
    }
    fprintf(bus_file, "%d", T_max);

    // Multiples Dutchführen der Simulation (minimaler Sinn)
    for (int i = 0; i < 3; i++) {
        // Kranksetzen der Busteilnehmer nach der Grundwahrscheinlichkeit sick_0
        for (int i_N = 0; i_N < N; i_N++) {
            if (gsl_rng_uniform(rng) <= sick_0) {
                bus_array[i_N] = 1;
            } else {
                bus_array[i_N] = 0;
            }

        }

        fprintf(bus_file, "\n%d", how_sick_is_the_bus(bus_array, N));
        printf("Schnupfen im Bus: %d\n", how_sick_is_the_bus(bus_array, N));

        // Itieren über T_max Zeitschritte
        for (int i_T = 0; i_T < T_max; i_T++) {
            int seat_n = gsl_rng_uniform(rng) * 40;                 // zufällig ausgewählter Sitz im Bus
            if (gsl_rng_uniform(rng) <= p) {                        // Spontangenesung der Wahrscheinlichkeit p
                bus_array[seat_n] = 0;
            } else { 
                if (bus_array[seat_n] == 1) {                       // Ansteckung der Nachbarn                    
                    if (seat_n == 0) {                              // erster Platz
                        bus_array[1] = 1;
                    } else if (seat_n == N-1) {                     // letzter Platz
                        bus_array[N-2] = 1;
                    } else {                                        // alle anderen Plätze
                        bus_array[seat_n - 1] = 1;
                        bus_array[seat_n + 1] = 1;
                    }
                } else {
                    printf("DICK\n");
                }
            }
            fprintf(bus_file, ", %d", how_sick_is_the_bus(bus_array, N));
            printf("Schnupfen im Bus: %d\n", how_sick_is_the_bus(bus_array, N));
        }
    }
    fclose(bus_file);
    return 0;
}