#ifndef NUMERICS_H
#define NUMERICS_H


#define cvc_PI 3.14159265358979323846


// Struktur eines 2er Tupels
struct cvc_tupel_2 {
    double x1;
    double x2;
};


// Typ gewöhnliche DGL
typedef int cvc_ode_func(double, const double[], double[], void*);
typedef int cvc_sde_func(double, const double[], double[], double, void*);


// Potenz für natürliche Zahlen x^n
double cvc_npow(double x, int n);

// Norm in 2 / 3 dimensions
double cvc_norm_2D(double x, double y);
double cvc_norm_3D(double x, double y, double z);


// Numerical Integration
double cvc_integrate(double x, int N);
double cvc_integrate_trapez(double left, double right, int N, double integrand(double));
double cvc_integrate_simpson(double func(double), double x, int N);
double cvc_integrate_simpson_2_param(double left, double right, double dx, double func(double, void*), void *params);

// Numerical Differentiation
double diff(double x, double delta, double func(double));

// Implementations of the error function 
double cvc_erf_simpson(double x, double delta_x);
double cvc_erf_midpoint(double x, double delta_x);


// Find root of given function func with respective parameters
double cvc_find_root_bisection(double func(double), double a, double b, double epsilon, int max_iter);
double cvc_find_root_newton_raphson(double func(double), double x0, double delta, double rel_tol, int max_iter);

// Solver of quadratic equations
struct cvc_tuple_2 cvc_solve_quadratic(double a, double b, double c);


// Tuple of 2 gaussian-distributed random numbers with static random-number-generator: Polarmethode
struct cvc_tuple_2 cvc_random_gaussian(void);

// MC-Berechnung der Dichte im Hyperquader [ai, bi] mit dim Dimensionen für die Dichtefunktion func()  
double cvc_mc_integrate(double func(double*, int), double a[], double b[], int dim, int N);

// 2-Dimensional Integration: Midpoint
double cvc_integrate_midpoint_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, double delta_x, double f(double, double));

// 2-Dimensional Integration: Monte-Carlo
double cvc_mc_integrate_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double));


// Numerical Integration using Euler / Runke-Kutta / Verlet methods with state array y and given parameters
void cvc_euler_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);
void cvc_rk2_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);
void cvc_rk4_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);
void cvc_verlet_step(double t, double delta_t, double y[], cvc_ode_func func, int dimension, void *params);

#endif