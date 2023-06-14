#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Parameter der Oszillation
const double OSZ_A0 = 2;
const double OSZ_T = 2000;
const double OSZ_OMEGA = 2 * cvc_PI / OSZ_T;


// Parameter des Ozeans
const double OCEAN_HEIGHT = 5000;                               // generelle Wassertiefe des Ozeans
const double OCEAN_HEIGHT_RIFF = 20;                            // Wassertiefe am Riff

const double OCEAN_X_LEFT = -1000e3;                            // linker Rand des Ozeans
const double OCEAN_X_RIGHT = 0;                                 // rechter Rand des Ozeans
const double OCEAN_X_RIFF_LEFT = -700e3;                        // linker Rand des Riffs
const double OCEAN_X_BEACH_LEFT = -550e3;                       // linker Rand des Strandes
const double OCEAN_X_BEACH_RIGHT = -450e3;                      // rechter Rand des Strandes
const double OCEAN_X_RIFF_RIGHT = -300e3;                       // rechter Rand des Riffs


// Integrationsparameter
const double T_MAX = 22000;
const double DELTA_X = 1e4;
const double DELTA_T = 10;
const int N = (OCEAN_X_RIGHT - OCEAN_X_LEFT) / DELTA_X + 1;


// Struktur der Parameter des d'Alembert-Integrators
struct params_dAlembert {
    double delta_t;
    double delta_x;
};


// Profil des Meeresbodes
double height_ocean(double x) {
    if (OCEAN_X_LEFT <= x && x < OCEAN_X_RIFF_LEFT) {
        return OCEAN_HEIGHT;
    } else if (OCEAN_X_RIFF_LEFT <= x && x < OCEAN_X_BEACH_LEFT) {
        return OCEAN_HEIGHT + (OCEAN_HEIGHT_RIFF-OCEAN_HEIGHT) * cvc_npow(sin(cvc_PI/2 * (x-OCEAN_X_RIFF_LEFT) / (OCEAN_X_BEACH_LEFT-OCEAN_X_RIFF_LEFT)), 2);
    } else if (OCEAN_X_BEACH_LEFT <= x && x < OCEAN_X_BEACH_RIGHT) {
        return OCEAN_HEIGHT_RIFF;
    } else if (OCEAN_X_BEACH_RIGHT <= x && x < OCEAN_X_RIFF_RIGHT) {
        return OCEAN_HEIGHT_RIFF + (OCEAN_HEIGHT_RIFF-OCEAN_HEIGHT) * cvc_npow(sin(cvc_PI/2 * (OCEAN_X_BEACH_RIGHT-x) / (OCEAN_X_BEACH_RIGHT-OCEAN_X_RIFF_RIGHT)), 2);
    } else if (OCEAN_X_RIFF_RIGHT <= x && x < OCEAN_X_RIGHT) {
        return OCEAN_HEIGHT;
    } else {
        return 0;
    }
}


// Anregung am freien Ende eines Seiles
double wave_excitation_A(double t) {
    return OSZ_A0 * sin(OSZ_OMEGA * t);
}


// Wellengeschwindigkeit für ein gegebenes x
double wave_speed_u(double x) {
    // return sqrt(cvc_EARTH_GRAVITATION * height_ocean(x));
    return 1;
}


// Berechnet den Zustand der Welle des nächsten Zeitschrittees
int dAlembertFTCS(double t, double y[], struct params_dAlembert) {
    double delta_x = params_dAlembert.delta_x;
    double delta_t = params_dAlembert.delta_t;

    // Berechnung der Wellengeschwindigkeit für x, x_(i-1) und x_(i+1)
    double u_i = wave_speed_u(x);
    double u_i_minus = wave_speed_u(x - delta_x);
    double u_i_plus = wave_speed_u(x + delta_x);

    // Berechnung des Parameters beta^2 aus den Wellengeschwindigkeiten für x, x_(i-1) und x_(i+1)
    double beta_i_squared = cvc_npow(u_i * delta_t / delta_x, 2);
    double beta_i_minus_squared = cvc_npow(u_i_minus * delta_t / delta_x, 2);
    double beta_i_plus_squared = cvc_npow(u_i_plus * delta_t / delta_x, 2);

    // Berechnung des Zustandes für den nächsten Zeitschritt
    double *y_temp = (double*) malloc(N * sizeof(double));
    for(size_t i = 0; i < N; i++) {
        double x = OCEAN_X_LEFT + (i * delta_x);

        double f = y[i % N];
        double f_west = y[(i - 1) % N];
        double f_east = y[(i + 1) % N];

        // if (0 <= t && t <= 2*cvc_PI/OSZ_OMEGA) {
        //     double f_west_ghost = A(t);                                                     // linker Ghost
        // } else {
        //     double f_west_ghost = 0;
        // }

        double f_west_ghost = A(t);                                                         // linker Ghost
        double f_east_ghost = 0;                                                            // rechter Ghost
        double f_t_minus = f + 1.0/2 * beta_i_squared * (f_west - 2*f + f_east);            // Berechnung f_(t-1) [sehr dubios!!!!!]

        double delta_f_t = cvc_diff

        if (i == 0) {
            y_temp[i] = 2*(1-beta_squared)*f - f_t_minus + beta_i_squared*(f_east_ghost+f_west) + 1.0/4*(beta_i_plus_squared-beta_i_minus_squared)*(f_east_ghost-f_west);
        } else if (i == N - 1) {
            y_temp[i] = 2*(1-beta_squared)*f - f_t_minus + beta_i_squared*(f_east+f_west_ghost) + 1.0/4*(beta_i_plus_squared-beta_i_minus_squared)*(f_east-f_west_ghost);
        } else {
            y_temp[i] = 2*(1-beta_squared)*f - f_t_minus + beta_i_squared*(f_east+f_west) + 1.0/4*(beta_i_plus_squared-beta_i_minus_squared)*(f_east-f_west);
        }
    }

    // Übertragen der neuen Werte in das Ursprungsaray
    for(int i = 0; i < N; i++) {
        y[i] = y_temp[i];
    }

    free(y_temp);
    return 0;
}


int main(void) {

    // Initialisieren des Ozeans
    double *y_wave = (double*) calloc(N, sizeof(double));
    double x = OCEAN_X_LEFT;                                                                // Laufvariable der Position im Ozean
    double t = 0;                                            

    FILE* ocean_file = fopen("data/A26_Tsunami_data.csv", "w");
    fprintf(ocean_file, "%g", x);
    for (int n_x = 0; n_x < N) {
        fprintf(ocean_file, ", %g", (n_x+1) * DELTA_X);
    }

    // Übergabestruktur von den d'Alembert-Integrator
    struct params_dAlembert = params_struct;

    while (t < T_MAX) {
        params_struct.delta_t = DELTA_T;
        params_struct.delta_x = DELTA_X;

        dAlembertFTCS(t, y_wave, params_struct);  


        t += DELTA_T;
    }
    return 0;
}