#include "definitions.h"
#include "support.h"

/*****************************************************************************************
MPI implementation of 1D Godunov scheme for fluid dynamics
******************************************************************************************/


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *R,				// density
		*P,					// pressure
		*U,					// velocity
		*RU,				// moment of impulse
		*RE,				// total energy
		*S;					// entropy

	double *FR,				// density flux
		*FRU,				// moment flux
		*FRE,				// energy flux
		*UFLUX;				// velocity flux

	double *uss, *pss, *dss;			// gas values on the boundaries
	double *x_init, *x_n, *x_n1;		// coordinates
	double *u_max_array;

	int iter = 0, count = 0;
	bool last = false;

	double timer, time_max, tau, dx, dtdx, len, x, CFL;
	double u1 = 0, u2 = 0, u3 = 0, u_loc = 0, u_max = 0;
	double start_t = 0, end_t = 0, duration = 0;
	double loop_time = 0, bound_time = 0, cfl_time = 0;

	/* For HALO MPI */
	double left_cell[3], right_cell[3];

	MPI_Status statuses[3];
	MPI_Request requests[3];


	/*** Start ***/

	/* Set number of cells */
	int numcells = 300;

	if (rank == 0)
	{
		printf("\nNumcells %d\n", numcells); fflush(0);
	}

	/* Information */
	if (rank == 0)
	{
		printf("List of nodes:\n"); fflush(0);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* Print processor name */
	char node_name[MPI_MAX_PROCESSOR_NAME];
	int node_name_length;
	MPI_Get_processor_name(node_name, &node_name_length);
	printf("%s\n", node_name); fflush(0);

	/* Read omp threads from the environment */
	int OMP_CORES = omp_get_max_threads();

	/* Chunk for OMP threading */
	int omp_chunk = numcells / OMP_CORES;

	/* Domain */
	len = LENGTH;
	dx = len / double(numcells);	// step 

	/* Create arrays */
	int loc_num = numcells / size;
	int rem = numcells % size;
	int stride1 = loc_num;
	int stride2 = 0;
	int r1 = rank;
	int r2 = 0;

	/* Set correct stride for dividing the computational
	domain if there is a reminder */
	if (rem != 0 && rank >= size / 2)
	{
		loc_num += 1;
		r1 = size / 2;
		r2 = rank - r1;
		stride1 = loc_num - 1;
		stride2 = loc_num;
	}
	
	mem_alloc(loc_num, &R, 32);
	mem_alloc(loc_num, &P, 32);
	mem_alloc(loc_num, &U, 32);
	mem_alloc(loc_num, &S, 32);
	mem_alloc(loc_num, &RU, 32);
	mem_alloc(loc_num, &RE, 32);

	mem_alloc(loc_num, &x_n, 32);
	mem_alloc(loc_num, &x_n1, 32);
	mem_alloc(loc_num, &x_init, 32);

	mem_alloc(loc_num + 1, &FR, 32);
	mem_alloc(loc_num + 1, &FRU, 32);
	mem_alloc(loc_num + 1, &FRE, 32);
	mem_alloc(loc_num + 1, &UFLUX, 32);
	mem_alloc(loc_num + 1, &dss, 32);
	mem_alloc(loc_num + 1, &uss, 32);
	mem_alloc(loc_num + 1, &pss, 32);

	u_max_array = (double*)malloc(size * sizeof(double));

	/*********************************************************************/

	start_t = MPI_Wtime();

	/* Mesh */
#pragma omp for simd schedule(simd:static)
		for (int i = 0; i < loc_num; i++)
		{
			x_init[i] = (r1 * stride1 + r2 * stride2 + i + 0.5)*dx;          // that are middles of cells
			x_n[i] = x_init[i];
		}

	/* Initial conditions */
	if (rank == 0) printf("Loop 0\n"); fflush(0);
#pragma omp parallel
	{ 
#pragma omp for simd schedule(simd:static)
		for (int i = 0; i < loc_num; i++)
		{
			R[i] = initial_density(x_init[i]);
			P[i] = initial_pressure(x_init[i]);
			U[i] = initial_velocity(x_init[i]);
		}
	}

	/* Computation of RU and RE */
	if (rank == 0) printf("Loop 1\n"); fflush(0);
#pragma omp parallel
	{
#pragma omp for simd schedule(simd:static)
		for (int i = 0; i < loc_num; i++)
		{
			RU[i] = R[i] * U[i];
			RE[i] = P[i] / (GAMMA - 1.0) + 0.5 * R[i] * U[i] * U[i]; // full mechanic energy
		}
	}

	/* Main loop of computational algrithm */
	iter = 0;
	timer = 0.0;
	time_max = time_max_array[PROBLEM];

	while (timer < time_max)
	{
		iter++;
		
#pragma omp parallel firstprivate(u1,u2,u3,u_loc) shared(u_max)
		{
			/* CFL condition */
#pragma omp for schedule(static, omp_chunk) nowait
			for (int i = 0; i < loc_num; i++)
			{
				u1 = U[i] + sqrt(GAMMA*P[i] / R[i]);
				u2 = U[i] - sqrt(GAMMA*P[i] / R[i]);
				u3 = U[i];

				if (u1 > u_loc) u_loc = u1;
				if (u2 > u_loc) u_loc = u2;
				if (u3 > u_loc) u_loc = u3;
			}
#pragma omp critical
			if (u_loc > u_max) u_max = u_loc;
		}

		/* Gather all u_max at rank 0 */
		MPI_Gather(&u_max, 1, MPI_DOUBLE, u_max_array, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/* Handling u_max at rank 0 */
		if (rank == 0)
		{
			for (int i = 0; i < size; i++)
				u_max = max(u_max_array[i], u_max);
		}

		/* Send u_max to the same variable at other processes */
		MPI_Bcast(&u_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
		

		/* Set time step via CFL number */
		tau = CFL04*dx / u_max;

		/* Handling the last time step */
		if (timer + tau > time_max)
		{
			tau = time_max - timer; // if the last time's step is bigger than distance to "time_max"
			last = true;	
		}
		dtdx = tau / dx;

		/* Boundary conditions */
		if (rank == 0) linear(R[0], U[0], P[0], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
		if (rank == size - 1) linear(R[loc_num - 1], U[loc_num - 1], P[loc_num - 1], R[loc_num - 1], U[loc_num - 1], P[loc_num - 1], dss[loc_num], uss[loc_num], pss[loc_num]);


		/* Nonlinear or linear Godunov scheme */
		if (last && rank == 0) printf("Loop 2\n"); fflush(0);
		if (timer < CROSS_POINT)
		{
			nonlinear_solver(loc_num, R, U, P, dss, uss, pss);
		}
		else
		{
			linear_solver(loc_num, R, U, P, dss, uss, pss);
		}

		/* HALO */

		int right = rank + 1;
		int left = rank - 1;

		// Messages should have the same tag!!!
		if (rank < size - 1)
		{
			MPI_Send(&R[loc_num - 1], 1, MPI_DOUBLE, right, 11, MPI_COMM_WORLD);
			MPI_Send(&U[loc_num - 1], 1, MPI_DOUBLE, right, 12, MPI_COMM_WORLD);
			MPI_Send(&P[loc_num - 1], 1, MPI_DOUBLE, right, 13, MPI_COMM_WORLD);

			MPI_Recv(&right_cell[0], 1, MPI_DOUBLE, right, 21, MPI_COMM_WORLD, &statuses[0]);
			MPI_Recv(&right_cell[1], 1, MPI_DOUBLE, right, 22, MPI_COMM_WORLD, &statuses[1]);
			MPI_Recv(&right_cell[2], 1, MPI_DOUBLE, right, 23, MPI_COMM_WORLD, &statuses[2]);

			linear(R[loc_num - 1], U[loc_num - 1], P[loc_num - 1], right_cell[0], right_cell[1], right_cell[2], dss[loc_num], uss[loc_num], pss[loc_num]);
		}
		
		if (rank > 0)
		{
			MPI_Recv(&left_cell[0], 1, MPI_DOUBLE, left, 11, MPI_COMM_WORLD, &statuses[0]);
			MPI_Recv(&left_cell[1], 1, MPI_DOUBLE, left, 12, MPI_COMM_WORLD, &statuses[1]);
			MPI_Recv(&left_cell[2], 1, MPI_DOUBLE, left, 13, MPI_COMM_WORLD, &statuses[2]);

			MPI_Send(&R[0], 1, MPI_DOUBLE, left, 21, MPI_COMM_WORLD);
			MPI_Send(&U[0], 1, MPI_DOUBLE, left, 22, MPI_COMM_WORLD);
			MPI_Send(&P[0], 1, MPI_DOUBLE, left, 23, MPI_COMM_WORLD);
		
			linear(left_cell[0], left_cell[1], left_cell[2], R[0], U[0], P[0], dss[0], uss[0], pss[0]);	
		}
		
		/* Computation of flux variables */
		if (last && rank == 0) printf("Loop 3\n"); fflush(0);
#pragma omp parallel
		{
#pragma omp for simd schedule(simd:static) 
			for (int i = 0; i <= loc_num; i++)
			{
				UFLUX[i] = uss[i];
				FR[i] = dss[i] * uss[i];
				FRU[i] = dss[i] * uss[i] * uss[i] + pss[i];
				FRE[i] = (pss[i] / (GAMMA - 1.0) + 0.5*dss[i] * uss[i] * uss[i])*uss[i] + pss[i] * uss[i];
			}
		}

		/* Euler coordinates */
#pragma omp parallel for simd schedule(simd:static)
		for (int i = 0; i < loc_num; i++)
		{
			x_n1[i] = x_n[i] + tau * UFLUX[i];
		}

		/* Computation of conservations laws (in conservative variables!) */
		if (last && rank == 0) printf("Loop 4\n"); fflush(0);
#pragma omp parallel
		{
#pragma omp for simd schedule(simd:static)
			for (int i = 0; i < loc_num; i++)
			{
				R[i] = R[i] - dtdx * (FR[i + 1] - FR[i]);
				RU[i] = RU[i] - dtdx * (FRU[i + 1] - FRU[i]);
				RE[i] = RE[i] - dtdx * (FRE[i + 1] - FRE[i]);
			}
		}

		/* Over-computation of velocity, pressure and entropy */
		if (last && rank == 0) printf("Loop 5\n"); fflush(0);
#pragma omp parallel
		{
#pragma omp for simd schedule(simd:static)
			for (int i = 0; i < loc_num; i++)
			{
				U[i] = RU[i] / R[i];
				P[i] = (GAMMA - 1.0) * (RE[i] - 0.5 * RU[i] * U[i]);
				S[i] = log(P[i] / pow(R[i], GAMMA));
			}
		}

		/* Addition to timer */
		timer += tau;

		/* Euler coordinates */
#pragma omp parallel for simd schedule(simd:static)
		for (int i = 0; i < loc_num; i++)
		{
			x_n[i] = x_n1[i];
		}

		MPI_Barrier(MPI_COMM_WORLD);

	} /******************************************* The end of iteration**************************/

	end_t = MPI_Wtime();
	duration = end_t - start_t;

	  /* Output to several files */
	output_last_step(loc_num, x_n, rank, R, U, P);

	mem_free(&R);
	mem_free(&P);
	mem_free(&U);
	mem_free(&S);
	mem_free(&RU);
	mem_free(&RE);

	mem_free(&FR);
	mem_free(&FRU);
	mem_free(&FRE);
	mem_free(&UFLUX);
	mem_free(&dss);
	mem_free(&uss);
	mem_free(&pss);

	mem_free(&x_n1);
	mem_free(&x_n);
	mem_free(&x_init);

	if (rank == 0) printf("\nSUCCESS!\nElapsed time: %lf sec\n", duration);
	
	MPI_Finalize();
	return 0;
}
