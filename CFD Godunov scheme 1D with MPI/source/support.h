#pragma once
/*************************
Header file for templates
**************************/


// Kernel
void sample(double &pm, double &um, double &s, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr, double &d, double &u, double &p);
void prefun(double &f, double &fd, double &p, double &dk, double &pk, double &ck);
double guessp(double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr);
void starpu(double &p, double &u, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr);

void nonlinear_solver(int numcells, double *R, double *U, double *P, double *dss, double *uss, double *pss);
void linear_solver(int numcells, double *R, double *U, double *P, double *dss, double *uss, double *pss, double** LOOP_TIME, int last);
void flux_count(FILE* *array_flux, int iter, int numcells, double timer, double tau, double *t, double *UFLUX);
void boundary_conditions(int numcells, double *dss, double *uss, double *pss, double *R, double *U, double *P);
void mem_alloc(int numcells, double **arr, int align);
void mem_free(double **arr);

double initial_density(double x);
double initial_pressure(double x);
double initial_velocity(double x);

double gyugonio(double p1, double ro1, double p2);
double sw_speed2(double ro1, double u1, double p1, double ro2, double p2);
double sw_speed(double ro1, double ro2, double u1, double u2);
double* finite_difference(int numb, double *mas);
double RW_prop(int digit, double x, double numb, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r);

void iteration(int numb, double* F_ro, double* ITER_TIME);
void linear(double dl, double ul, double pl, double dr, double ur, double pr, double &d, double &u, double &p);
void linear_check(double dl, double ul, double pl, double dr, double ur, double pr, int &left, int &middle, int &right, int numb);
void runge(double *massiv, int ldf, int numb);
void rw_diff_num_analit(int numb, int numcells, double *R, double *U, double *P);
void analitical_RW(FILE* file_name, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r, double numb);
void analitical_SW(int numcells, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r, double *res_p, double *res_u, double *res_d, double timer);
void analitical_riemann(int numcells, double p1, double ro1, double u1, double p2, double ro2, double u2, double *sol_p, double *sol_u);
void analitical_riemann_modeling(int numcells, double ro1, double u1, double p1, double ro2, double u2, double p2, double timer, double *all_d, double *all_u, double *all_p);
void analitical_writing_into_file(int numcells, double *R_D, double *R_U, double *R_P, double timer);
void difference_analitical_riemann_Linf(int numb, double *R, double *U, double *P, double *R_D, double *R_U, double *R_P, double &delta_ro, double &delta_u, double &delta_p);
void difference_analitical_riemann_L1(int numb, double *R, double *U, double *P, double *R_D, double *R_U, double *R_P, double &sum_ro, double &sum_u, double &sum_p);
void difference_SW(int numcells, double timer, double *R, double *U, double *P, double *shw_diff_d, double *shw_diff_u, double *shw_diff_p, double *shw_analit_d, double *shw_analit_u, double *shw_analit_p);
void outline_integral_riemann(int numcells, double timer, double tau, const double tt1, const double tt2, double xx1, double xx2, double *xx, double* R, double*U, double*P, double*RE, double*S, /*output*/ double sum[4][4]);
void file_exact_diff(int numcells, double *exact_R, double *exact_U, double *exact_P, double *exact_RE, double *exact_S, double *diff_R, double *diff_U, double *diff_P, double time);
void inf_before_start(int numcells, double *R, double *U, double *P, double &D_analit);
void first_step_validation(FILE *file, int numcells, int c_c, double timer, double *R, double *U, double *P, double *dss, double *uss, double *pss);
void file_n_smooth_steps(int numcells, double timer, double tau, double *x_layer, double *R, double *U, double* P, double *RE, double *S, double *S_diff, double *UFLUX);
void null_array(double *arr, int a, int b);
void null_array(int *arr, int a, int b);
void set_bound(FILE *plot, int i);
void output_last_step(int numcells, double dx, double D_analit, double *R, double *U, double *P);


// Gnuplot
void gnuplot_one_iteration(int numcells);
void gnuplot_RW_DIFF(int numcells);
void gnuplot_RW_NUM_ANALITIC(int numcells);
void gnuplot_P_PLUS_PG(int numcells);
void gnuplot_all_iter_one_time(int run[NUM_ITER], int numb, double time);
void gnuplot_all_iterations_NC(int numb);
void gnuplot_one_it_NC();
void gnuplot_conservative(int numb);
void gnuplot_five_t_steps(int numb);
void gnuplot_n_smooth_NC(int numb);
void gnuplot_n_smooth_NC2(int numcells, int *n_r, int *n_u, int *n_p);
void gnuplot_analitical_riemann(int numcells, double *R, double *U, double *P, double *R_D, double *R_U, double *R_P);
void gnuplot_n_smooth(int numb);
void gnuplot_n_smooth2(int numcells, int sw1[3][N_smooth], int sw2[3][N_smooth], int sw3[3][N_smooth]);
void gnuplot_n_smooth3(int numcells);
void gnuplot_analitical_riemann2(int numcells, int *n_r, int *n_u, int *n_p);
