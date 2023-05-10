#ifndef NUMERICS_H
#define NUMERICS_H


#define CVC_PI 3.14159265358979323846

// Typ gewöhnliche DGL
typedef int ode_func(double, const double[], double[], void*);

// Potenz für natürliche Zahlen x^n
double npow(double x, int n);

// Norm in 2 dimensions
double norm_2D(double x, double y);


// Numerical Integration
double integrate(double x, int N);
double integrate_trapez(double left, double right, int N, double integrand(double));
double integrate_simpson(double func(double), double x, int N);
double integrate_simpson_2_param(double left, double right, double dx, double func(double, void*), void *params);

// Numerical Differentiation
double diff(double x, double delta, double func(double));

// Implementations of the error function 
double erf_simpson(double x, double delta_x);
double erf_midpoint(double x, double delta_x);


// Find root of given function func with respective parameters
double find_root_bisection(double func(double), double a, double b, double epsilon, int max_iter);
double find_root_newton_raphson(double func(double), double x0, double delta, double rel_tol, int max_iter);

// Solver of quadratic equations
struct tuple_2 solve_quadratic(double a, double b, double c);


// Tuple of 2 gaussian-distributed random numbers with static random-number-generator: Polarmethode
struct tuple_2 random_gaussian(void);

// MC-Berechnung der Dichte im Hyperquader [ai, bi] mit dim Dimensionen für die Dichtefunktion func()  
double mc_integrate(double func(double*, int), double a[], double b[], int dim, int N);

// 2-Dimensional Integration: Midpoint
double integrate_midpoint_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, double delta_x, double f(double, double));

// 2-Dimensional Integration: Monte-Carlo
double mc_integrate_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double));


// Numerical Integration using Euler / Runke-Kutta / Verlet methods with state array y and given parameters
void euler_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);
void rk2_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);
void rk4_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);
void verlet_step(double t, double delta_t, double y[], ode_func func, int dimension, void *params);

#endif