#pragma once
/*****************************
Header file for preprocessor 
definitions and extern pattern 
of arrays from arrays.cpp

Each definition correspond to some
additional job performed during
work of the scheme.

User can uncomment needed #define
to enable the respective job.
*****************************/
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <omp.h>
#include <mpi.h>
#include <xmmintrin.h>

/************** MACRO ****************/

// 1) General preprocessor definition

//#define PRINT printf("Printing with GNUPLOT is set up\n");
//#define OUTPUT_N_SMOOTH printf("Output to several files. One file = one layer")
//#define INTEGRAL printf("Computation of accuracy via using counter integrals")

#if (PROBLEM < 3 || PROBLEM == 7 || PROBLEM == 4)
//#define FLUX_COUNT
#endif

#ifdef INTEGRAL
#define RUNGE
#undef PRINT
#endif

// 2) Other processing

//#define NC printf("Coordinates of shockwave is setted up\n");  // It influences on output in .dat file
//#define NC2 // (x-x0-Dt)/h
//#define P_PLUS_PG printf("Printing new task about pressure and pressure gradient");
//#define SW_FINITE_DIFF  printf("Its a work with shock wave and right parts. We are determinating the wavelength\n")
//#define RW_NUM_ANALITICAL  printf("We are conculating the difference between numeric and analitical solvers\n") 
//#define FIVE_T_STEPS
//#define DIFF_ANALIT_RIEMANN
//#define BOOST
//#define L1_NORM
//#define ENTROPY_CHECK
//#define SW_POINTS_PRINT
//#define CFL_SWITCH
//#define DIFF_ANALYT
//#define BOUND_COND
//#define DEBUG


#define NUM_ITER 8

#define GRID 3

#define N_smooth 40
#define N_bound 50

#define LOOPS 6         // the number of loops under the time calculations

#define EPS		1.0e-6	// the threshold for the solving nonlinear solution
#define MAX_ITER	20	// the number of iteration to find nonlinear solution

#define PROBLEM		2
/*
0 - shock wave
1 - rarify wave
2 - shock tube
3 - periodic continious
4 - shock tube with 2 shock waves
5 - two shock tube (распад разрыва) в x=4 и x=6
6 - rarifaction wave in x=6 comes to shockwave and intersects it in x=7;
7 - two rarefaction waves propagates in different directions
8 - the same as problem 6 + count of points on the waves; strong sw
9 - test Sod
10 - some experiments
11 - shock wave intersect with rarefaction domain
12 - test Toro 1
13 - test Toro 2
14 - test Toro 3
15 - test Toro 4
16 - test Toro 5
17 - shock tube [0:10], x = 5 - discontinuity
18 - shock wave with boundary condition in x = 0
19 - shock wave with boundary condition U=u(t)
20 - shock wave intersect another shock wave. Bound condition
*/

/*************************************/

/* Accuracy */
#ifdef INTEGRAL
#define X1 0.05
#define X2 0.95
#define T1 0.0
#define T2 0.5
#endif


/* Switch between nonlinear and linear cases */

//#define EXACT_DISC
#ifdef EXACT_DISC
#define CROSS_POINT time_max_array[PROBLEM]
#define TYPE 'E'
#else
//#define CROSS_POINT 1.875  problem 1
#define CROSS_POINT 0.000
#define TYPE 'L'
#endif

/*************************************/

#if (PROBLEM == 8 )
#define DELTA 0.005
#else
#define DELTA 0.005
#endif

#define CFL04 0.4
#define CFL08 0.8

#define PI			3.1415926535897932

#define A_TERM		1.0			// term = - A * \vec{u} h^{2K-1} (dp/dx)^{2K}
#define K_TERM		2.0			

#define GAMMA 1.4
#define R0 8.93
#define C0 3.97

#define C1 1.0


/**************************************/

#define X_G(x) (pow((x), GAMMA))
#define SQ_2(x) ((x)*(x)/2)

/**************************************/

#define RUNGE_KUTTA	0	/*	1 - Runge-Kutta method
0 - Classic method
*/

/***************************************/
#if (PROBLEM==4)
#define P4_ONE_WAVE
#endif

// SW 1 // SW 2 // SW 3

#define st_R3 1.0
#define st_P3 1.0
#define st_U3 0.0


#define st_P2 1.4017
#define st_R2 gyugonio(st_P3, st_R3, st_P2)
#define st_U2 ((sw_speed2(st_R3, st_U3, st_P3, st_R2, st_P2)*(1.0 - (st_R3)/(st_R2))) + ((st_R3)*(st_U3)/(st_R2)))

#ifdef P4_ONE_WAVE
#define st_P1 5.0
#else
#define st_P1 3.288303 //5
#endif 
#define st_R1 gyugonio(st_P2, st_R2, st_P1)
#define st_U1 ((sw_speed2(st_R2, st_U2, st_P2, st_R1, st_P1)*(1.0 - (st_R2)/(st_R1))) + ((st_R2)*(st_U2)/(st_R1)))

/**************************************************/

#define st_th_P2 1.0
#define st_th_R2 1.0
#define st_th_U2 0.0

#define st_th_P1 3.5
#define st_th_R1 gyugonio(st_th_P2, st_th_R2, st_th_P1)
#define st_th_U1 ((sw_speed2(st_th_R2, st_th_U2, st_th_P2, st_th_R1, st_th_P1)*(1.0 - (st_th_R2)/(st_th_R1))) + ((st_th_R2)*(st_th_U2)/(st_th_R1)))

/* Define point of discontinuity */
#if (PROBLEM==0 || PROBLEM==1)
#define DISC_POINT 0.1         // 0.3 - rarification wave test
#elif (PROBLEM == 2 || PROBLEM == 9)
#define DISC_POINT 0.5
#elif (PROBLEM == 4)
#define DISC_POINT 0.5
#elif (PROBLEM == 5)
#define DISC_POINT 0
#elif (PROBLEM == 6)
#define DISC_POINT 1.0
#elif (PROBLEM == 7)
#define DISC_POINT 0.5
#elif (PROBLEM == 8) || (PROBLEM==20)
#define DISC_POINT 0.2
#elif (PROBLEM == 10)
#define DISC_POINT 0.2
#elif (PROBLEM == 11)
#define DISC_POINT 0.2
#elif (PROBLEM == 12)
#define DISC_POINT 0.3
#elif (PROBLEM == 13)
#define DISC_POINT 0.5
#elif (PROBLEM == 14)
#define DISC_POINT 0.5
#elif (PROBLEM == 15)
#define DISC_POINT 0.4
#elif (PROBLEM == 16)
#define DISC_POINT 0.8
#elif (PROBLEM == 17)
#define DISC_POINT 5.0
#elif (PROBLEM == 18)
#define DISC_POINT 0.0
#elif (PROBLEM == 19)
#define DISC_POINT 0.0
#else
#define DISC_POINT 0.0
#endif

/* Define length of a problem */
#if (PROBLEM == 5 || PROBLEM == 17)
#define LENGTH 10.0
#elif (PROBLEM == 19)
#define LENGTH 20.0
#else
#define LENGTH 1.0
#endif

extern char prop[];
extern char dip[];
extern double time_max_array[];
extern float percents[NUM_ITER];

#define max(a,b) ((a) > (b)) ? (a) : (b)

// arrays.cpp

extern int nmesh[];
extern int nprnt[];

extern const double g1, g2, g3, g4, g5, g6, g7, g8;

extern float left_SW[];
extern float right_SW[];

extern float left_19[];
extern float right_19[];

extern float left_7[];
extern float right_7[];

extern float left_SW_cs[];
extern float right_SW_cs[];

extern float left_ST[];
extern float right_ST[];

extern float left_RW[];
extern float right_RW[];

extern float left_sodd[];
extern float right_sodd[];

extern float left_2RR[];
extern float right_2RR[];

extern float left_RW_SW[];
extern float right_RW_SW[];

extern float left_RW_RW[];
extern float right_RW_RW[];

extern float left_RW1_RW1[];
extern float right_RW1_RW1[];

extern float left_RW2_RW2[];
extern float right_RW2_RW2[];

extern float left_RW3_RW3[];
extern float right_RW3_RW3[];

extern float left_RP[];
extern float right_RP[];

extern float left_SW_SW[];
extern float right_SW_SW[];

// iteration.cpp
// R U P

extern double delta_RUP[3][NUM_ITER];
extern double l1_RUP[3][NUM_ITER];



