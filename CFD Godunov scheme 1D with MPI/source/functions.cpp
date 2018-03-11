#include "definitions.h"
#include "support.h"

/****************************************
Source file contains linear and
nonlinear solvers of the Riemann problem
and additional functionalities for
processing results of the computations
****************************************/

/* SIMD linear solver */
void linear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss, double** LOOP_TIME, int last)
{
	double wtime = 0;
	double *C = new double[numcells];
	double *RC = new double[numcells];
	double *H = new double[numcells];

#pragma omp parallel private(wtime)
	{
		wtime = omp_get_wtime();

#pragma omp for simd schedule(simd:static) // schedule(dynamic,64) 
		for (int i = 1; i < numcells; i++)
		{

			C[i - 1] = sqrt(GAMMA*P[i - 1] / R[i - 1]);
			C[i] = sqrt(GAMMA*P[i] / R[i]);

			RC[i - 1] = R[i - 1] * C[i - 1];
			RC[i] = R[i] * C[i];

			H[i - 1] = 1.0 / (RC[i - 1]);
			H[i] = 1.0 / (RC[i]);

			if (U[i - 1] > C[i - 1])
			{
				pss[i] = P[i - 1];
				uss[i] = U[i - 1];
				dss[i] = R[i - 1];

			}
			else if (U[i] < -C[i])
			{
				pss[i] = P[i];
				uss[i] = U[i];
				dss[i] = R[i];
			}
			else
			{
				pss[i] = (U[i - 1] - U[i] + P[i - 1] * H[i - 1] + P[i] * H[i]) / (H[i - 1] + H[i]);
				uss[i] = (RC[i - 1] * U[i - 1] + RC[i] * U[i] + P[i - 1] - P[i]) / (RC[i - 1] + RC[i]);

			    if (uss[i] > 0) dss[i] = R[i - 1] - R[i - 1] / C[i - 1] * (uss[i] - U[i - 1]);
			    else dss[i] = R[i] + R[i] / C[i] * (uss[i] - U[i]);

			}

		}
		wtime = omp_get_wtime() - wtime;
		LOOP_TIME[2][omp_get_thread_num()] += wtime;
		if (last) printf("Time taken by thread %d is %f\n", omp_get_thread_num(), LOOP_TIME[2][omp_get_thread_num()]);
	}

	delete[] RC;
	delete[] C;
	delete[] H;
}

/* Scalar linear solver */
void linear(double dl, double ul, double pl, double dr, double ur, double pr, double &d, double &u, double &p) // передача параметров по ссылке
{
	double bigU, bigP, bigS, bigR, cl, cr, help, hl, hr, R3, R4;

	cl = sqrt(GAMMA*pl / dl);
	cr = sqrt(GAMMA*pr / dr);
	hl = 1.0 / (dl*cl);
	hr = 1.0 / (dr*cr);

	if (ul > cl)
	{
		bigP = pl;
		bigU = ul;
		bigR = dl;

	}
	else if (ur < -cr)
	{
		bigP = pr;
		bigU = ur;
		bigR = dr;
	}
	else
	{
		bigP = (ul - ur + pl / (dl*cl) + pr / (dr*cr)) / (hl + hr);
		bigU = (dl*cl*ul + dr*cr*ur + pl - pr) / (dl*cl + dr*cr);
		/*	if (bigU >= 0) bigS = pl / pow(dl, GAMMA);
		else bigS = pr / pow(dr, GAMMA);
		help = bigP / bigS;
		bigR = pow(help, 1.0 / GAMMA);*/
		R3 = dl - dl / cl * (bigU - ul);
		R4 = dr + dr / cr * (bigU - ur);
		if (bigU > 0) bigR = R3;
		else bigR = R4;
	}

	u = bigU;
	p = bigP;
	d = bigR;

}


void nonlinear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss)
{

	double um = 0, pm = 0;
	double s_char;

	/***** Параметры ударной трубы для соседних ячеек *****/
	// Параметры слева
	double	dl,  // плотность
		ul,  // скорость
		pl,  // давление
		cl;  // скорость звука

			 // Параметры справа
	double	dr,  // плотность
		ur,  // скорость
		pr,  // давление
		cr;  // скорость звука

	//====================================== EXACT RIEMANN PROBLEM =========================
	for (int i = 1; i < numcells; i++)
	{
		// левая сторона
		dl = R[i - 1];
		pl = P[i - 1];
		ul = U[i - 1];
		cl = sqrt(GAMMA*pl / dl);
		// правая сторона
		dr = R[i];
		pr = P[i];
		ur = U[i];
		cr = sqrt(GAMMA*pr / dr);

		starpu(pm, um, dl, ul, pl, cl, dr, ur, pr, cr);

		// Решение задачи распада разрыва
		sample(pm, um, s_char, dl, ul, pl, cl, dr, ur, pr, cr, dss[i], uss[i], pss[i]);
	}
}

void boundary_conditions(int numcells, double *dss, double *uss, double *pss, double *R, double *U, double *P)
{
	double c0, s0, cn, sn;

#if (PROBLEM==18 || PROBLEM==20)
	// set pressure
	pss[0] = initial_pressure(0.0);
	//	pss[numcells] = initial_pressure(LENGTH);

	c0 = sqrt(GAMMA*P[0] / R[0]);
	//	cn = sqrt(GAMMA*P[numcells] / R[numcells]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	//double rn_const = U[numcells] + P[numcells] / (R[numcells] * cn);

	uss[0] = l0_const + pss[0] / (R[0] * c0);
	//	uss[numcells] = rn_const - pss[numcells] / (R[numcells] * cn);

	s0 = log(P[0] / pow(R[0], GAMMA));
	//	sn = log(P[numcells] / pow(R[numcells], GAMMA));

	//dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);
	//	dss[numcells] = pow(pss[numcells] / sn, 1.0 / GAMMA);

	//R3 = dl - dl / cl * (bigU - ul);
	dss[0] = R[0] + R[0] / c0 * (uss[0] - U[0]);

	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);

#elif (PROBLEM == 19)

	uss[0] = timer;
	c0 = sqrt(GAMMA*P[0] / R[0]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	double rn_const = U[numcells] + P[numcells] / (R[numcells] * Cn);

	pss[0] = (uss[0] - l0_const)*(R[0] * c0);

	s0 = log(P[0] / pow(R[0], GAMMA));
	dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);


#elif (PROBLEM == 3)
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[numcells], uss[numcells], pss[numcells]);
#else
	linear(R[0], U[0], P[0], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);
#endif
}

void flux_count(FILE* *array_flux, int iter, int numcells, double timer, double tau, double *t, double *UFLUX)
{
	int t_ind[N_bound] = { 0 };
	int numcells_flux;
	numcells_flux = numcells;

	double dx = LENGTH / double(numcells);

	for (int i = 0; i < N_bound; i++)
	{
		t_ind[i] = i * numcells_flux / N_bound;
	}

	if (iter == 1)
	{
		for (int i = 0; i < N_bound; i++)
			t[i] = (t_ind[i] + 0.5)*dx;
	}
	else
	{
		for (int i = 0; i < N_bound; i++)
		{
			//t[i] = t[i] + UFLUX[t_ind[i]] * tau;

			fprintf(array_flux[i], "%lf %lf %lf\n", t[i], timer, UFLUX[t_ind[i]]);
		}
	}
		//t[i] = (t_ind[i] + 0.5)*dx - UFLUX[t_ind[i]] * timer;
}

/**************************************************
Получение плотности, давления и скорости в точке
**************************************************/
void sample(double &pm, double &um, double &s, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr,
	double &d, double &u, double &p)
{

	double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

	if (s <= um)
	{
		// точка слева от контактного разрыва
		if (pm <= pl)
		{
			// левая волна разрежения
			shl = ul - cl;

			if (s <= shl)
			{
				// точка слева
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				cml = cl * pow(pm / pl, g1);
				stl = um - cml;

				if (s > stl)
				{
					// точка слева от контактного разрыва
					d = dl * pow(pm / pl, 1.0 / GAMMA);
					u = um;
					p = pm;
				}
				else
				{
					// точка в волне разрежения
					u = g5 * (cl + g7*ul + s);
					c = g5 * (cl + g7*(ul - s));
					d = dl*pow(c / cl, g4);
					p = pl*pow(c / cl, g3);
				}
			}
		}
		else
		{
			// левая ударная волна
			pml = pm / pl;
			sl = ul - cl * sqrt(g2*pml + g1);

			if (s <= sl)
			{
				// точка слева от ударной волны
				d = dl;
				u = ul;
				p = pl;
			}
			else
			{
				// точка слева от контактного разрыва
				d = dl * (pml + g6) / (pml*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
	}
	else
	{
		// точка справа от контактного разрыва
		if (pm > pr)
		{
			// правая ударная волна
			pmr = pm / pr;
			sr = ur + cr * sqrt(g2*pmr + g1);

			if (s >= sr)
			{
				// точка справа от ударной волны
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				// точка справа от контактного разрыва
				d = dr * (pmr + g6) / (pmr*g6 + 1.0);
				u = um;
				p = pm;
			}
		}
		else
		{
			// правая волна разрежения
			shr = ur + cr;

			if (s >= shr)
			{
				// точка справа от волны разрежения
				d = dr;
				u = ur;
				p = pr;
			}
			else
			{
				cmr = cr * pow(pm / pr, g1);
				str = um + cmr;

				if (s <= str)
				{
					// точка справа от контактного разрыва
					d = dr * pow(pm / pr, 1.0 / GAMMA);
					u = um;
					p = pm;
				}
				else
				{
					// точка в левой волне разрежения
					u = g5 * (-cr + g7*ur + s);
					c = g5 * (cr - g7*(ur - s));
					d = dr * pow(c / cr, g4);
					p = pr * pow(c / cr, g3);
				}
			}
		}
	}

}

/**************************************************
Решение для одной из частей
**************************************************/
void prefun(double &f, double &fd, double &p,
	double &dk, double &pk, double &ck)
{
	double ak, bk, pratio, qrt;

	if (p <= pk)
	{
		// волна разрежения
		pratio = p / pk;
		f = g4*ck*(pow(pratio, g1) - 1.0);
		fd = (1.0 / (dk*ck))*pow(pratio, -g2);

	}
	else
	{
		// ударная волна
		ak = g5 / dk;
		bk = g6*pk;
		qrt = sqrt(ak / (bk + p));
		f = (p - pk)*qrt;
		fd = (1.0 - 0.5*(p - pk) / (bk + p)) * qrt;
	}

}


/**************************************************
Вычисление начального приближения
**************************************************/
double guessp(double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr)
{
	double	cup, gel, ger,
		pmax, pmin, ppv, pq, pm,
		ptl, ptr,
		qmax, quser, um;

	/*** Вычисление приближения давления из PVRS решения Римана ***/
	quser = 2.0;
	cup = 0.25 * (dl + dr) * (cl + cr);
	ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
	ppv = ppv > 0.0 ? ppv : 0.0;
	pmin = pl > pr ? pr : pl;
	pmax = pl > pr ? pl : pr;
	qmax = pmax / pmin;

	if (qmax <= quser && pmin <= ppv && ppv <= pmax)
	{
		pm = ppv;
	}
	else
	{
		if (ppv < pmin)
		{
			// две волны разрежения
			pq = pow(pl / pr, g1);
			um = (pq*ul / cl + ur / cr + g4*(pq - 1.0)) / (pq / cl + 1.0 / cr);
			ptl = 1.0 + g7*(ul - um) / cl;
			ptr = 1.0 + g7*(um - ur) / cr;
			pm = 0.5*(pl*pow(ptl, g3) + pr*pow(ptr, g3));
		}
		else
		{
			// две ударных волны
			gel = sqrt((g5 / dl) / (g6*pl + ppv));
			ger = sqrt((g5 / dr) / (g6*pr + ppv));
			pm = (gel*pl + ger*pr - (ur - ul)) / (gel + ger);
		}

	}
	return pm;
}


/**************************************************
Определение давления и скорости контактного разрыва
**************************************************/
void starpu(double &p, double &u, double dl, double ul, double pl, double cl, double dr, double ur, double pr, double cr)
{
	int i;

	double pstart,	// начальное приближение давления 
		pold,	// предыдущее приближение давления
		udiff,	// разность скоростей
		change;	// разность давлений

	double	fl, fld, // функции давления
		fr, frd;

	/*** Вычисление начального приближения ***/
	pstart = guessp(dl, ul, pl, cl, dr, ur, pr, cr);

	/*** Предыдущее приближение ***/
	pold = pstart;

	/*** Разность скоростей ***/
	udiff = ur - ul;

	/*** Метод Ньютона для определения давления ***/

	// разность между разными приближениями давлений
	change = 10.0 * EPS;

	for (i = 0; i<MAX_ITER && change>EPS; i++)
	{
		// решение для левой части
		prefun(fl, fld, pold, dl, pl, cl);

		// решение для правой части
		prefun(fr, frd, pold, dr, pr, cr);

		// очередное приближение давления
		p = pold - (fl + fr + udiff) / (fld + frd);

		// разность между разными приближениями давлений
		change = 2.0 * fabs((p - pold) / (p + pold));

		// если давление отрицательное, до обнуляем его
		if (p < 0.0) p = 0.0;

		pold = p;
	}

	// определение скорости
	u = 0.5 * (ul + ur + fr - fl);
}

void mem_alloc(int numcells, double* *arr, int align)
{
	*arr = (double*)_mm_malloc(numcells * sizeof(double), align);
}

void mem_free(double **arr)
{
	_mm_free(*arr);
}

void first_step_validation(FILE *file3, int numcells, int c_c, double timer, double *R, double *U, double *P, double *dss, double *uss, double *pss)
{
	double len = LENGTH;
	double dx = len / double(numcells);
	double x;


	fprintf(file3, "cc: %d timer: %lf\n\n", c_c, timer);
	
	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;
		fprintf(file3, "%lf %lf %lf %lf %lf %lf %lf\n", x, R[i], U[i], P[i], dss[i], uss[i], pss[i]);
	}
	fprintf(file3,"*********************************************************\n\n");
	
}


void linear_check(double dl, double ul, double pl, double dr, double ur, double pr, int &left, int &middle, int &right, int numb)
{
	double cl, cr, help, hl, hr;
	cl = sqrt(GAMMA*pl / dl);
	cr = sqrt(GAMMA*pr / dr);
	hl = 1.0 / (dl*cl);
	hr = 1.0 / (dr*cr);

	if (numb == 1)
	{
		if (ul > cl || ur < -cr)
		{
			if (ul > cl)
			{
				left++;
			}
			if (ur < -cr)
			{
				right++;
			}
		}
		else
		{
			middle++;
		}
	}

	if (numb == 2)
	{
		if (ul + cl < 0 || ur - cr > 0)
		{
			if (ul + cl < 0)
			{
				left++;
			}
			if (ur - cr > 0)
			{
				right++;
			}
		}
		else
		{
			middle++;
		}
	}

}

double gyugonio(double p1, double ro1, double p2/*за ударной волной*/)
{
	return ro1*((GAMMA + 1.0)*p2 + (GAMMA - 1.0)*p1) / ((GAMMA - 1.0)*p2 + (GAMMA + 1.0)*p1);
	//  R2 = ro2*((GAMMA + 1.0)*P + (GAMMA - 1.0)*p2) / ((GAMMA - 1.0)*P + (GAMMA + 1.0)*p2);
}

double sw_speed(double ro1, double ro2, double u1, double u2)
{
	return (ro1*u1 - ro2*u2) / (ro1 - ro2);

}

double sw_speed2(double ro1, double u1, double p1, double ro2 /*за ударной волной*/, double p2 /*за ударной волной*/)
{
	return u1 + sqrt(ro2*(p2 - p1) / ro1 / (ro2 - ro1));
}

double law(double t)
{
	return sqrt(t);
}


double initial_density(double x)
{
	switch (PROBLEM)
	{
	case 19:
	case 18:
	case 0:	if (x <= DISC_POINT) return 1.271413930046081;
			else return 1.0;
	case 1:	if (x <= DISC_POINT) return 1.551608179649565;
			else return 2.0;
	/*case 1:	if (x <= DISC_POINT)
		return 1.0;
			else
				return 0.585507;*/
	case 2:	if (x <= 0.5)
		return 2.0;
			else
				return 1.0;
	case 3:	return 1.0;
		/*	case 4: if (x <= DISC_POINT) return 1.271413930046081;
		else if (x >= 1.0 - DISC_POINT) return 1.271413930046081;
		else return 1.0;*/
	case 4: if (x <= DISC_POINT) return 3.0;
			else return 2.0;

	case 5: if (x <= 3.0) return 6;
			else if (x >= 7.0) return 6;
			else return 4.0;
	case 6: if (x <= 0.2) return 1.36205;
			else if (x <= 0.3) return 4.875;
			else return 3;
	case 7: if (x < DISC_POINT) return 1.0;
			else return 1.0; // 2.0
	case 8: if (x <= 0.05) return st_R1;
			else if (x <= 0.3) return st_R2;
			else return st_R3;
	case 9: if (x <= DISC_POINT) return 1.0;
			else return 0.125;
	case 10: if (x <= DISC_POINT) return st_th_R1;
			 else return st_th_R2;
			 //case 10: return -0.55*tanh(2 * PI*(x - 0.5)) + 1.546; // weak
			 //case 10: return -0.9*tanh(2 * PI*(x - 0.5)) + 1.9; // strong
	case 11: if (x <= 0.1) return st_R2;
			 else if (x <= 0.6) return st_R3;
			 else 0.2;
		// TORO tests
		// 1
	case 12: if (x <= DISC_POINT) return 1.0;
			 else return 0.125;
		// 2
	case 13: return 1.0;
		// 3
	case 14: return 1.0;
		// 4
	case 15: if (x <= DISC_POINT) return 5.99924;
			 else return 5.99242;
		// 5
	case 16: return 1.0;
	case 17:	if (x <= 5.0) return 2.0;
				else return 1.0;
	case 20: if (x <= 0.00) return st_R1;
			else if (x <= DISC_POINT) return st_R2;
			else return st_R3;
	}
	return 0;
}

double initial_pressure(double x)
{
	switch (PROBLEM)
	{
	case 19:
	case 18:
	case 0:	if (x <= DISC_POINT) return 1.401789770179879;
			else return 1.0;
	case 1:	if (x <= DISC_POINT) return 1.401789770179879;
			else return 2.0;
	/*case 1:	if (x <= DISC_POINT)
		return 1.0;
			else
				return 0.466727;*/

	case 2:	if (x <= DISC_POINT) return 2.0;
			else return 1.0;
	case 3:	return 1.0;
		/*case 4: if (x <= DISC_POINT) return 1.401789770179879;
		else if (x >= 1.0 - DISC_POINT) return 1.401789770179879;
		else return 1.0;*/
	case 4: if (x <= DISC_POINT) return 2.0;
			else return 1.0;

	case 5: if (x <= 3.0) return 7;
			else if (x >= 7.0) return 7;
			else return 2.0;
	case 6: if (x <= 0.2) return 0.67099;
			else if (x <= 0.3) return 4;
			else return 2;
	case 7: return 1.0;
	case 8: if (x <= 0.05) return st_P1;
			else if (x <= 0.3) return st_P2;
			else return st_P3;
	case 9: if (x <= DISC_POINT) return 1.0;
			else return 0.1;
			//case 10: return -tanh(2 * PI*(x - 0.5)) + 2; //weak

			//case 10: return -2 * tanh(2 * PI*(x - 0.5)) + 3; //strong


	case 10: if (x <= DISC_POINT) return st_th_P1;
			 else return st_th_P2;
	case 11: if (x <= 0.1) return st_P2;
			 else return st_P3;
		// TORO tests
		// 1
	case 12: if (x <= DISC_POINT) return 1.0;
			 else return 0.1;
		// 2
	case 13: return 0.4;
		// 3
	case 14: if (x <= DISC_POINT) return 1000.0;
			 else return 0.01;
		// 4
	case 15: if (x <= DISC_POINT) return 460.894;
			 else return 46.095;
		// 5
	case 16: if (x <= DISC_POINT) return 1000.0;
			 else return 0.01;
	case 17: if (x <= 5.0) return 2.0;
			 else return 1.0;
	case 20: if (x <= 0.00) return st_P1;
			else if (x <= DISC_POINT) return st_P2;
			else return st_P3;
	}
	return 0;
}

double initial_velocity(double x)
{
	switch (PROBLEM)
	{
	case 19:
	case 18:
	case 0:	if (x <= DISC_POINT) return 0.292868067614595;
			else return 0.0;

	case 1:	if (x <= DISC_POINT) return -0.292868067614595;
			else  return 0.0;

	/*case 1:	if (x <= DISC_POINT)
		return 0.75 + 0.4533;// -0.43;
			else
				return 1.386 + 0.4533;// -0.43;*/

	case 2:	return 0.0;
		//if (x <= 0.5)  2 ударных волны
		//	return 1.0;
		//else return 0.0;
	case 3:	return sin(2.0*PI*x);
		/*case 4: if (x <= DISC_POINT) return 0.292868067614595;
		else if (x >= 1.0 - DISC_POINT) return 0.292868067614595;
		else return 0.0;*/
	case 4: if (x <= DISC_POINT) return 4.0;
			else return 2.0;

	case 5: if (x <= 3.0) return 0;
			else if (x >= 7.0) return 0.0;
			else return 0.0;
	case 6: if (x <= 0.2) return -0.7;
			else if (x <= 0.3) return 0.50637;
			else return 0;
	case 7: if (x <= DISC_POINT) return -1.0;
			else return 1.0;
	case 8: if (x <= 0.05) return st_U1;
			else if (x <= 0.3) return st_U2;
			else return st_U3;
	case 9: return 0.0;
	case 10: if (x <= DISC_POINT) return st_th_U1;
			 else return st_th_U2;
			 //	case 10: return -0.51*tanh(2 * PI*(x - 0.5)) + 0.51; // weak

			 //case 10: return -0.8*tanh(2 * PI*(x - 0.5)) + 0.8; // strong
	case 11: if (x <= 0.1) return st_U2;
			 else return st_U3;
		// TORO tests
		// 1
	case 12: if (x <= DISC_POINT) return 0.75;// +0.2;// +0.4533;
			 else return 0.0;// +0.2;// 0.4533;
		// 2
	case 13: if (x <= DISC_POINT) return -2.0;
			 else return 2.0;
		// 3
	case 14: return 0.0;
		// 4
	case 15: if (x <= DISC_POINT) return 19.5975;
			 else return -6.19633;
		// 5
	case 16: return -19.59745;
	case 17: return 0.0;
	case 20: if (x <= 0.00) return st_U1;
			else if (x <= DISC_POINT) return st_U2;
			else return st_U3;
	}
	return 0;
}

void file_n_smooth_steps(int numcells, double timer, double tau, double *x_layer,
	double *R, double *U, double* P, double *RE, double *S, double *S_diff, double *UFLUX)
{
	FILE* fout, *fout_NC;
	char FileName[255], FileName2[255];
	double dx = LENGTH / double(numcells);
	double ps, us, ds, cs, es, es_d;
	double D_analit = 1.371913;

	int proverka[N_smooth] = { 0 };
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;

	for (int i = 0; i < N_smooth; i++)
	{
		time_control[i] = (i + 1)*k_step;
	}

	for (int k = 0; k < N_smooth; k++)
	{

		if (fabs(timer - time_control[k]) <= tau && proverka[k] == 0)
		{
#ifndef NC					
			sprintf(FileName, "workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4lf.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (char)TYPE, time_control[k]);
			fout = fopen(FileName, "w");
			fprintf(fout, "Timer: %lf\n", timer);
#else
			sprintf(FileName2, "workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4lf_NC.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (char)TYPE, time_control[k]);
			fout_NC = fopen(FileName2, "w");
			fprintf(fout_NC, "Timer: %lf\n", timer);
#endif
			for (int i = 0; i < numcells; i++)  // вывод всегда по 100 точек ( с первой итерации которые )
			{

				ds = R[i];
				//us = U[i] - 0.4533;
				us = U[i];
				ps = P[i];
				cs = sqrt(GAMMA*P[i] / R[i]);
				es = S[i];
				es_d = S_diff[i];

#ifndef NC
				fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf %lf\n", x_layer[i], ds, us, ps, cs, es, es_d);
#else
				fprintf(fout_NC, "%9.6lf %lf %lf %lf %lf %lf %lf\n", x, ds, us, ps, cs, es, es_d);
#endif
			}
#ifndef NC
			fclose(fout);
#else
			fclose(fout_NC);
#endif
			proverka[k] = 1;

		}
	}
}


void difference_SW(int numcells, double timer, double *R, double *U, double *P, 
	                          double *shw_diff_d, double *shw_diff_u, double *shw_diff_p, 
	                          double *shw_analit_d, double *shw_analit_u, double *shw_analit_p)
{
	/*******************difference analit and numeric solutions************/
	/**********************shock wave*************************************/

	analitical_SW(numcells, initial_pressure(0.05), initial_density(0.05), initial_velocity(0.05),
		                    initial_pressure(0.2), initial_density(0.2), initial_velocity(0.2),
		                    shw_analit_p, shw_analit_u, shw_analit_d, timer);
	for (int j = 0; j < numcells; j++)
	{
	shw_diff_p[j] = P[j] - shw_analit_p[j];
	shw_diff_u[j] = U[j] - shw_analit_u[j];
	shw_diff_d[j] = R[j] - shw_analit_d[j];
	}
}

void analitical_RW(FILE* file_name, double ip_l, double id_l, double iu_l,
	double ip_r, double id_r, double iu_r, double numb)
{
	FILE* out_ul, *out_pl, *out_dl;
	FILE* out_ur, *out_pr, *out_dr;
	out_pl = fopen("data_pl.txt", "w");
	out_ul = fopen("data_ul.txt", "w");
	out_dl = fopen("data_dl.txt", "w");

	out_pr = fopen("data_pr.txt", "w");
	out_ur = fopen("data_ur.txt", "w");
	out_dr = fopen("data_dr.txt", "w");

	double xl, xr, c, l0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	//tt = 1 / (u2 + c2 - D);
	tt = numb*0.1;   // only this one have NUMB count since 1  !!!

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	fprintf(file_name, "u(x)=%lf*(x-0.1)/%lf+(%lf)\n", q3, tt, q11);

	xl = (iu_l - q11)*tt / q3 + 0.1;
	xr = (iu_r - q11)*tt / q3 + 0.1;
	fprintf(file_name, "xu_l=%lf\nxu_r=%lf\n", xl, xr);
	fprintf(out_ul, "0 %lf\n%lf %lf", iu_l, xl, iu_l);
	fprintf(out_ur, "%lf %lf\n1 %lf", xr, iu_r, iu_r);

	// для давления  p=(q4*x-q5)^q6

	q4 = q1 / (sqrt(GAMMA)*pow(A, 1 / (2 * GAMMA)));
	q5 = q4*l0;
	q6 = q2*GAMMA;
	fprintf(file_name, "p(x)=(%lf*(x-0.1)/%lf-(%lf))**%lf\n", q4, tt, q5, q6);  // почему, если разрвы в 0.1 мы должны отнять

	xl = (pow(ip_l, (1 / q6)) + q5)*tt / q4 + 0.1;
	xr = (pow(ip_r, (1 / q6)) + q5)*tt / q4 + 0.1;
	fprintf(file_name, "xp_l=%lf\nxp_r=%lf\n", xl, xr);

	fprintf(out_pl, "0 %lf\n%lf %lf", ip_l, xl, ip_l);
	fprintf(out_pr, "%lf %lf\n1 %lf", xr, ip_r, ip_r);

	// для плотности ro=(q7*x-q8)^q2

	q7 = q1 / sqrt(GAMMA*A);
	q8 = q7*l0;
	fprintf(file_name, "ro(x)=(%lf*(x-0.1)/%lf-(%lf))**%lf\n", q7, tt, q8, q2);

	xl = (pow(id_l, (1 / q2)) + q8)*tt / q7 + 0.1;
	xr = (pow(id_r, (1 / q2)) + q8)*tt / q7 + 0.1;
	fprintf(file_name, "xd_l=%lf\nxd_r=%lf\n\n", xl, xr);

	fprintf(out_dl, "0 %lf\n%lf %lf", id_l, xl, id_l);
	fprintf(out_dr, "%lf %lf\n1 %lf", xr, id_r, id_r);

	fprintf(file_name, "id_l(x)=%lf\n", id_l);
	fprintf(file_name, "ip_l(x)=%lf\n", ip_l);
	fprintf(file_name, "iu_l(x)=%lf\n\n", iu_l);

	fprintf(file_name, "id_r(x)=%lf\n", id_r);
	fprintf(file_name, "ip_r(x)=%lf\n", ip_r);
	fprintf(file_name, "iu_r(x)=%lf\n", iu_r);

	fclose(out_dl);
	fclose(out_dr);
	fclose(out_ul);
	fclose(out_ur);
	fclose(out_pl);
	fclose(out_pr);

	return;
}

void analitical_SW(int numcells, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r, double *res_p, double *res_u, double* res_d, double timer)
{
	double D_analit;
	double x;
	double *x_lay = new double[numcells];
	double dx = LENGTH / double(numcells);
	for (int i = 0; i < numcells; i++)
		x_lay[i] = 0;

	D_analit = (id_r * iu_r - id_l * iu_l) / (id_r - id_l);

	x = D_analit * timer + 0.1; // discontinuity point

	for (int i = 0; i < numcells; i++)
	{
		x_lay[i] = i*dx + 0.5*dx;
		if (x_lay[i] <= x)
		{
			res_p[i] = ip_l;
			res_d[i] = id_l;
			res_u[i] = iu_l;
		}
		else
		{
			res_p[i] = ip_r;
			res_d[i] = id_r;
			res_u[i] = iu_r;
		}
	}

	free(x_lay);

}

/*modeling of exact solution of riemann problem*/
void analitical_riemann(int numcells, double p1, double ro1, double u1, double p2, double ro2, double u2, double *sol_p, double *sol_u)
{
	//double P = 0, U = 0;  //искомые давление и скорость (постоянные значения слева и справа от контактного разрыва)
	double R1, E1; //искомые значения плотности и внутренней энергии СЛЕВА от контактного разрыва
	double R2, E2; //искомые значения плотности и внутренней энергии СПРАВА от контактного разрыва

	// начальные данные

	double c1; // давление, плотность, скорость СЛЕВА
	double c2; // давление, плотность, скорость СПРАВА

	c1 = sqrt(GAMMA*p1 / ro1);
	c2 = sqrt(GAMMA*p2 / ro2);

	double a1; // массовая скорость слева
	double a2; // массовая скорость справа

	double P_prev = 0;
	double P_now = 0;

	double delta = 0;
	double drob = 0;
	double h1 = 0, h2 = 0;
	int l = 0;  //счетчик

	double fi_P1 = 0;
	double alfa = 0;
	double help = 0, h3 = 0, h4 = 0;
	double z = 0;

	P_prev = 20; // P_prev > p1 и P_prev > p2;    ударные волны СЛЕВА и СПРАВА

	do
	{
		if (P_prev > p1) a1 = sqrt(ro1*((GAMMA + 1)*P_prev / 2.0 + (GAMMA - 1)*p1 / 2.0));
		else
		{
			h2 = (GAMMA - 1.0) / (2.0 * GAMMA);
			h1 = P_prev / p1;
			drob = (1.0 - P_prev / p1) / (1.0 - pow(h1, h2));
			a1 = (GAMMA - 1.0) / (2.0 * GAMMA)*ro1*c1*drob;
		}
		drob = 0;

		if (P_prev > p2) a2 = sqrt(ro2*((GAMMA + 1)*P_prev / 2.0 + (GAMMA - 1)*p2 / 2.0));
		else
		{
			h2 = (GAMMA - 1) / (2.0 * GAMMA);
			h1 = P_prev / p2;
			drob = (1.0 - P_prev / p2) / (1.0 - pow(h1, h2));
			a2 = (GAMMA - 1) / (2.0 * GAMMA)*ro2*c2*drob;
		}
		drob = 0;

		z = P_prev / (p1 + p2);
		h3 = (GAMMA + 1) / 2.0 / GAMMA;
		h4 = (GAMMA - 1) / 2.0 / GAMMA;
		help = ((GAMMA - 1) / 3.0 / GAMMA)*((1.0 - z) / (pow(z, h3)*(1.0 - pow(z, h4)))) - 1;
		if (help > 0) alfa = help;
		else alfa = 0;
		fi_P1 = (a2*p1 + a1*p2 + a1*a2*(u1 - u2)) / (a1 + a2);
		P_now = (alfa*P_prev + fi_P1) / (1 + alfa);

		delta = P_now - P_prev;

		// переприсваивание

		P_prev = P_now;
		P_now = 0;

		l++;

	} while (fabs(delta) > 0.0000001);

	*sol_p = P_prev;
	*sol_u = (a1*u1 + a2*u2 + p1 - p2) / (a1 + a2);
}

void rw_diff_num_analit(int numb, int numcells, double *R, double *U, double *P)
{
	/* Rarify wave - numeric vs analitic */
	double c, l0, A, tt;
	double iu_l, id_l, ip_l, iu_r, id_r, ip_r;
	double x, x_NC, xl, xr;
	double q1, q11, q2, q3;
	int check1 = 0, check2 = 0;

	double dx = LENGTH / double(numcells);

	iu_l = initial_velocity(0.05);
	id_l = initial_density(0.05);
	ip_l = initial_pressure(0.05);

	iu_r = initial_velocity(0.2);
	id_r = initial_density(0.2);
	ip_r = initial_pressure(0.2);

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	tt = (numb + 1)*0.1;

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	xl = (iu_l - q11)*tt / q3 + DISC_POINT; // счет по старому, по итерациям - какая итерация, такое и время. 
	xr = (iu_r - q11)*tt / q3 + DISC_POINT; // в нашем случае ЭТО НЕВЕРНО, так как теперь итерации отвечают только за ШАГ СЕТКИ, время должно быть ФИКСИРОВАНО
									 // здесь ошибка при вычислении точного решения!!!


									 // 0-U, 1-P, 2-R, 3-RU, 4-RE
	double** difference_RW;
	difference_RW = new double*[5];

	for (int j = 0; j < 5; j++)
		difference_RW[j] = new double[numcells];

	double* RW_R, *RW_P, *RW_U;
	RW_R = new double[numcells];
	RW_P = new double[numcells];
	RW_U = new double[numcells];

	int counter = 0, counter_all = 0, counter2 = 0;

	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;

		if (x < xl)
		{
			RW_U[i] = iu_l; // FOR RARIFY WAVE
			RW_P[i] = ip_l;
			RW_R[i] = id_l;
			difference_RW[0][i] = U[i] - RW_U[i];
			difference_RW[1][i] = P[i] - RW_P[i];
			difference_RW[2][i] = R[i] - RW_R[i];
			difference_RW[3][i] = 0;
			difference_RW[4][i] = 0;
			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter2++;
			}

		}
		if (x >= xl && x <= xr)
		{
			counter_all++;
			RW_U[i] = RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[0][i] = U[i] - RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			RW_P[i] = RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[1][i] = P[i] - RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter++;
				counter2++;
			}

			RW_R[i] = RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[2][i] = R[i] - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);

			difference_RW[3][i] = R[i] * U[i] - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)*RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r);
			difference_RW[4][i] = R[i] * (P[i] / pow(R[i], GAMMA) + SQ_2(U[i])) - RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)*(RW_prop(1, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r) / pow(RW_prop(2, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r), GAMMA) + SQ_2(RW_prop(0, x, numb, ip_l, id_l, iu_l, ip_r, id_r, iu_r)));
		}
		if (x > xr)
		{
			RW_U[i] = iu_r;
			RW_P[i] = ip_r;
			RW_R[i] = id_r;

			difference_RW[0][i] = U[i] - RW_U[i];
			difference_RW[1][i] = P[i] - RW_P[i];
			difference_RW[2][i] = R[i] - RW_R[i];
			difference_RW[3][i] = 0;
			difference_RW[4][i] = 0;
			if (fabs(difference_RW[1][i]) <= 0.02)
			{
				counter2++;
			}

		}

	}

	check1 = 0;
	check2 = 0;

	int *i_helper = new int[10];
	double *x_helper = new double[10];

	double xl_num, xr_num;

	/****************Boundary of numerical rarify wave******************/

	int i_mem_left = 0, i_mem_right = 0;

	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;
		if (x >= xl && check1 == 0)
		{
			xl_num = x;   // по координате x
			i_mem_left = i;  // по счетчику i
			check1 = 1;
		}
		if (x >= xr && check2 == 0)
		{
			xr_num = x;
			i_mem_right = i;
			check2 = 1;
		}
	}
	printf("%lf %lf", xl_num, xr_num);
	/****************Boundary of numerical rarify wave******************/

	if (numb > 0)
	{
		int helper = counter_all / 10;

		for (int j = 0; j < 10; j++)
		{
			i_helper[j] = i_mem_left + j*helper;
		}
		printf("\n");
		for (int i = 0; i < numcells; i++)
		{
			x = i*dx + 0.5*dx;
			for (int j = 0; j < 10; j++)
			{
				if (i == i_helper[j])
				{
					x_helper[j] = x - (U[i] + sqrt(GAMMA*P[i] / R[i]))*time_max_array[PROBLEM];
					//		printf("%lf\n", x_helper[j]);
				}
			}
		}
	}

	printf("\nxl: %lf, xr: %lf\n", xl, xr);
	printf("points in rarify wave %d %d\n", counter, counter_all);
	printf("points2 %d %d\n", counter2, numcells);
	percents[numb] = float(counter) / float(counter_all) * 100.0f;  //when in int i counter is devided by counter all, its a деление нацело, so the result of 1/4=0;
	printf("percents of the middle: %f\n", percents[numb]);
	percents[numb] = float(counter2) / float(numcells) * 100.0f;  //when in int i counter is devided by counter all, its a деление нацело, so the result of 1/4=0;
	printf("percents of all stream [0:1]: %f\n", percents[numb]);


	FILE *out4;
	FILE *out5;
	char *FileName2, *FileName3;
	FileName2 = new char[64];
	FileName3 = new char[64];
	sprintf(FileName2, "N%04d_RW_difference.dat", numcells);
	sprintf(FileName3, "N%04d_RW_NUM_ANALITIC.dat", numcells);
	out4 = fopen(FileName2, "w");
	out5 = fopen(FileName3, "w");
	check1 = 0;
	check2 = 0;
	for (int i = 0; i < numcells; i++)
	{
#ifndef NC
		x = i*dx + 0.5*dx;
		if (x >= xl && check1 == 0)
		{
			fprintf(out5, "left b: %lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
			check1 = 1;
			continue;
		}
		if (x >= xr && check2 == 0)
		{
			fprintf(out5, "right b: %lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
			check2 = 1;
			continue;
		}
		fprintf(out4, "%lf %lf %lf %lf\n", x, difference_RW[0][i], difference_RW[1][i], difference_RW[2][i]);
		fprintf(out5, "%lf %lf %lf %lf %lf %lf %lf %lf\n", x, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i], U[i] + sqrt(GAMMA*P[i] / R[i]));
#else

		// Выведем в новых координатах
		x_NC = (i*dx + 0.5*dx - 0.1) / (time_max_array[0] * (numb + 1));
		fprintf(out4, "%lf %lf %lf %lf\n", x_NC, difference_RW[0][i], difference_RW[1][i], difference_RW[2][i]);
		fprintf(out5, "%lf %lf %lf %lf %lf %lf %lf\n", x_NC, U[i], RW_U[i], P[i], RW_P[i], R[i], RW_R[i]);
#endif
	}
	fclose(out4);
	fclose(out5);

}

void analitical_riemann_modeling(int numcells, double ro1, double u1, double p1, double ro2, double u2, double p2, double timer,
	/*output*/double *all_d, double *all_u, double *all_p)
{
	static double P, U; //solution of Riemann problem
	double R1, R2;
	int numcells2 = numcells / 2;
	double c1, C, prop, r0, A, tt;
	double q1, q2, q3, q4, q6, q7;
	int static ex_sol = 0;

	double v1, V;

	double *xx = new double[numcells];
	double dx = LENGTH / double(numcells);

	if (ex_sol == 0)
	{
		analitical_riemann(numcells, p1, ro1, u1, p2, ro2, u2, &P, &U); //для всех шагов сетки и шагов по времени нам нужно получить точное решение всего 1 раз
	}

	c1 = sqrt(GAMMA*p1 / ro1);
	r0 = u1 + 2 * c1 / (GAMMA - 1); // для левой волны разрежения; r0 инвариант - постоянен
	A = p1 / (pow(ro1, GAMMA));
	C = (r0 - U)*(GAMMA - 1) / 2.0;

	V = U - C;
	v1 = u1 - c1;
	//находим плотность R1 за волной разрежения через cвойство постоянности инварианта римана r0
	R1 = GAMMA*P / (C*C);
	//находим плотность sol_d[0] за ударной волной c помощью адиабаты Гюгонио
	R2 = ro2*((GAMMA + 1.0)*P + (GAMMA - 1.0)*p2) / ((GAMMA - 1.0)*P + (GAMMA + 1.0)*p2);

	/*контактный разрыв*/
	double x_KR;
	x_KR = U*timer + 0.5;

	/*Строим решение справа от x=0.5 - ударная волна*/
	double D_analit;
	double x;

	D_analit = (ro2* u2 - R2 * U) / (ro2 - R2);
	x = D_analit * timer + 0.5; // discontinuity point

	double c1_rw, C_rw, l0, A_rw;

	c1_rw = sqrt(GAMMA*p1 / ro1);
	l0 = u1 - 2 * c1_rw / (GAMMA - 1); // для правой волны разрежения; l0 инвариант - постоянен
	A_rw = p1 / (pow(ro1, GAMMA));
	C_rw = (U - l0)*(GAMMA - 1) / 2.0;

	if (ex_sol == 0)
	{
		printf("u1-c1: %lf\nU-C: %lf, D_analit: %lf, D_cont razr: %lf\n", v1, V, D_analit, U);
		printf("t=0.2: x1: %lf, x2: %lf, x3: %lf, x4: %lf\n", v1*0.2 + 0.5, V*0.2 + 0.5, D_analit*0.2 + 0.5, U*0.2 + 0.5);
		printf("rariry wave t=0.25: x1: %lf, x2: %lf\n", (u1 + c1_rw)*0.25 + DISC_POINT, (U + C)*0.25 + DISC_POINT);
		ex_sol = 1;
	}

	for (int i = 0; i < numcells; i++)
	{
		xx[i] = i*dx + 0.5*dx;
	}

	for (int i = numcells2; i < numcells; i++)
	{
		if (xx[i] <= x)
		{
			all_p[i] = P;
			all_u[i] = U;
		}
		else
		{
			all_p[i] = p2;
			all_d[i] = ro2;
			all_u[i] = u2;
		}
	}

	/*отдельно для плотности расчет, учитываем скорость контактного разрыва*/
	for (int i = 0; i < numcells; i++)
	{
		if (xx[i] >= x_KR && xx[i] <= x) all_d[i] = R2;
	}


	/*Строим решение слева от x=0.5 - волна разрежения*/
	double xl, xr;

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	q4 = q1 / (sqrt(GAMMA)*pow(A, 1 / (2 * GAMMA)));
	q6 = q2*GAMMA;
	q7 = q1 / sqrt(GAMMA*A);

	xl = (u1 - q1*r0)*timer / q3 + 0.5;
	xr = (U - q1*r0)*timer / q3 + 0.5;
	//	printf("%lf %lf\n", xl, xr);

	for (int i = 0; i < numcells2; i++)
	{
		if (xx[i] < xl)
		{
			all_p[i] = p1;
			all_d[i] = ro1;
			all_u[i] = u1;
		}
		else if (xx[i] >= xl && xx[i] <= xr)
		{
			all_p[i] = pow(q4*(r0 - (xx[i] - 0.5) / timer), q6);
			all_d[i] = pow(q7*(r0 - (xx[i] - 0.5) / timer), q2);
			all_u[i] = q3*(xx[i] - 0.5) / timer + q1*r0;
		}
		else
		{
			all_p[i] = P;
			all_u[i] = U;
		}
	}

	for (int i = 0; i < numcells; i++)
	{
		if (xx[i] >= xr && xx[i] <= x_KR) all_d[i] = R1;
	}

	free(xx);
}

void null_array(double *arr, int a, int b)
{
	for (int i = a; i < b; i++)
		arr[i] = 0.0;
}

void null_array(int *arr, int a, int b)
{
	for (int i = a; i < b; i++)
		arr[i] = 0;
}


void difference_analitical_riemann_Linf(int numb, double *R, double *U, double *P, double *R_D, double *R_U, double *R_P, double &delta_ro, double &delta_u, double &delta_p)
{
	double dx = LENGTH / double(nmesh[numb]);
	double *x_lay;
	double *difference_R;
	double *difference_U;
	double *difference_P;

	int new_numcells = int(0.41*nmesh[numb]);

	mem_alloc(nmesh[numb], &x_lay, 32);
	mem_alloc(new_numcells, &difference_R, 32);
	mem_alloc(new_numcells, &difference_U, 32);
	mem_alloc(new_numcells, &difference_P, 32);

	null_array(x_lay, 0, nmesh[numb]);
	delta_ro = 0;
	delta_u = 0;
	delta_p = 0;


#pragma omp parallel for simd schedule(simd: static)
	for (int i = 0; i < new_numcells; i++)
	{
		x_lay[i] = i*dx + 0.5*dx;
		difference_R[i] = fabs(R[i] - R_D[i]);
		difference_U[i] = fabs(U[i] - R_U[i]);
		difference_P[i] = fabs(P[i] - R_P[i]);
	}

	for (int i = 0; i < new_numcells; i++)
	{
		if (difference_R[i] > delta_ro) delta_ro = difference_R[i];
		if (difference_U[i] > delta_u) delta_u = difference_U[i];
		if (difference_P[i] > delta_p) delta_p = difference_P[i];
	}

	delta_ro = delta_ro * dx;
	delta_u = delta_u * dx;
	delta_p = delta_p * dx;

	mem_free(&x_lay);
	mem_free(&difference_P);
	mem_free(&difference_R);
	mem_free(&difference_U);
}

void difference_analitical_riemann_L1(int numb, double *R, double *U, double *P, double *R_D, double *R_U, double *R_P, double &sum_ro, double &sum_u, double &sum_p)
{
	double dx = LENGTH / double(nmesh[numb]);
	double *x_lay;
	double *difference_R;
	double *difference_U;
	double *difference_P;

	int new_numcells = int(0.41*nmesh[numb]);

	mem_alloc(nmesh[numb], &x_lay, 32);
	mem_alloc(new_numcells, &difference_R, 32);
	mem_alloc(new_numcells, &difference_U, 32);
	mem_alloc(new_numcells, &difference_P, 32);

#pragma omp parallel for simd schedule(simd: static)
	for (int i = 0; i < new_numcells; i++)
	{
		x_lay[i] = i*dx + 0.5*dx;
		difference_R[i] = fabs(R[i] - R_D[i]);
		difference_U[i] = fabs(U[i] - R_U[i]);
		difference_P[i] = fabs(P[i] - R_P[i]);
	}

#pragma omp parallel for reduction(+:sum_ro,sum_u,sum_p)
	for (int i = 0; i < new_numcells; i++)
	{
		sum_ro += difference_R[i];
		sum_u += difference_U[i];
		sum_p += difference_P[i];
	}

	sum_ro = sum_ro * dx;
	sum_u = sum_u * dx;
	sum_p = sum_p * dx;

	mem_free(&x_lay);
	mem_free(&difference_P);
	mem_free(&difference_R);
	mem_free(&difference_U);
}

/* 0 - U, 1 - P, 2 - D */
double RW_prop(int digit, double x, double numb, double ip_l, double id_l, double iu_l, double ip_r, double id_r, double iu_r)
{
	double c, prop, l0, A, tt;
	double q1, q2, q3, q4, q5, q6, q7, q8, q11;

	c = sqrt(GAMMA*ip_r / id_r);
	l0 = iu_r - 2 * c / (GAMMA - 1);
	A = ip_r / (pow(id_r, GAMMA));
	tt = (numb + 1)*0.1; 	//tt = 1 / (u2 + c2 - D);

	// для скорости u=q3*x+q11

	q1 = (GAMMA - 1) / (GAMMA + 1);
	q11 = q1*l0;
	q2 = 2 / (GAMMA - 1);
	q3 = 2 / (GAMMA + 1);
	if (digit == 0)
	{
		return prop = q3*(x - 0.1) / tt + q11;
	}

	// для давления  p=(q4*x-q5)^q6

	q4 = q1 / (sqrt(GAMMA)*pow(A, 1 / (2 * GAMMA)));
	q5 = q4*l0;
	q6 = q2*GAMMA;
	if (digit == 1)
	{
		return prop = pow((q4*(x - 0.1) / tt - q5), q6);

	}

	// для плотности ro=(q7*x-q8)^q2

	q7 = q1 / sqrt(GAMMA*A);
	q8 = q7*l0;

	if (digit == 2)
	{
		return prop = pow((q7*(x - 0.1) / tt - q8), q2);
	}

	return 0;
}

void outline_integral_riemann(int numcells, double timer, double tau, double tt1, double tt2, double xx1, double xx2, double *xx, double* R, double*U, double*P, double*RE, double*S,
	/*output*/ double sum[4][4])
{
	int static check1 = 0;
	int static check2 = 0;

	double static sum_l_M, sum_l_I, sum_l_S, sum_l_E;
	double static sum_r_M, sum_r_I, sum_r_S, sum_r_E;
	double static sum_t_M, sum_t_I, sum_t_S, sum_t_E;
	double static sum_b_M, sum_b_I, sum_b_S, sum_b_E;

	int static numcells_check;

	if (numcells == 100) numcells_check = 100;

	/*проверка на новый пространственный шаг*/
	if (numcells_check != numcells)
	{
		numcells_check = numcells;
		sum_l_M = 0, sum_l_I = 0, sum_l_S = 0, sum_l_E = 0;
		sum_r_M = 0, sum_r_I = 0, sum_r_S = 0, sum_r_E = 0;
		sum_t_M = 0, sum_t_I = 0, sum_t_S = 0, sum_t_E = 0;
		sum_b_M = 0, sum_b_I = 0, sum_b_S = 0, sum_b_E = 0;
		check1 = 0;
		check2 = 0;
	}

	int OMP_CORES = omp_get_max_threads();

	double dx = LENGTH / double(numcells);
	int omp_chuck = numcells / OMP_CORES;

	if (timer >= tt1 && timer <= tt2)
	{
		int l_bound = int((xx1 - 0.5*dx) / dx);
		int r_bound = int((xx2 - 0.5*dx) / dx);

		/*************massa*****************/

		sum_l_M += R[l_bound] * U[l_bound] * tau;
		sum_r_M += R[r_bound] * U[r_bound] * tau;

		/***********impulse****************/

		sum_l_I += (P[l_bound] + R[l_bound] * U[l_bound] * U[l_bound]) * tau;
		sum_r_I += (P[r_bound] + R[r_bound] * U[r_bound] * U[r_bound]) * tau;

		/***********entropy****************/

		sum_l_S += R[l_bound] * S[l_bound] * U[l_bound] * tau;
		sum_r_S += R[r_bound] * S[r_bound] * U[r_bound] * tau;

		/***********energy****************/

		sum_l_E += (RE[l_bound] + P[l_bound]) * U[l_bound] * tau;
		sum_r_E += (RE[r_bound] + P[r_bound]) * U[r_bound] * tau;
	}

	if (timer >= tt1 && check1 == 0)
	{

#pragma omp parallel for reduction(+:sum_b_M,sum_b_I,sum_b_S,sum_b_E) schedule(guided) num_threads(OMP_CORES)
		for (int i = 0; i < numcells; i++)
		{
			if (xx[i] >= xx1 && xx[i] <= xx2)
			{
				sum_b_M += R[i] * dx;
				sum_b_I += R[i] * U[i] * dx;
				sum_b_S += R[i] * S[i] * dx;
				sum_b_E += RE[i] * dx;
			}
		}

		check1 = 1;
	}


	if (timer >= tt2 && check2 == 0)
	{
#pragma omp parallel for reduction(+:sum_t_M,sum_t_I,sum_t_S,sum_t_E) schedule(guided) num_threads(OMP_CORES)
		for (int i = 0; i < numcells; i++)
		{
			if (xx[i] >= xx1 && xx[i] <= xx2)
			{
				sum_t_M += R[i] * dx;
				sum_t_I += R[i] * U[i] * dx;
				sum_t_S += R[i] * S[i] * dx;
				sum_t_E += RE[i] * dx;
			}
		}
		check2 = 1;
	}

	sum[0][0] = sum_t_M;
	sum[0][1] = sum_t_I;
	sum[0][2] = sum_t_S;
	sum[0][3] = sum_t_E;

	sum[1][0] = sum_b_M;
	sum[1][1] = sum_b_I;
	sum[1][2] = sum_b_S;
	sum[1][3] = sum_b_E;

	sum[2][0] = sum_r_M;
	sum[2][1] = sum_r_I;
	sum[2][2] = sum_r_S;
	sum[2][3] = sum_r_E;

	sum[3][0] = sum_l_M;
	sum[3][1] = sum_l_I;
	sum[3][2] = sum_l_S;
	sum[3][3] = sum_l_E;
}

void inf_before_start(int numcells, double *R, double *U, double *P, double &D_analit)
{ 
#if (PROBLEM == 0)
	D_analit = (R[numcells - 1] * U[numcells - 1] - R[0] * U[0]) / (R[numcells - 1] - R[0]);
	printf("Analitical speed: %10.8lf\n\n", D_analit);
#elif (PROBLEM == 12)
	double uc_left, uc_right;

	//	uc_left = U[0] +0.43 - sqrt(GAMMA*P[0] / R[0]);
	//uc_right = U[numcells] + 0.43  - sqrt(GAMMA*P[numcells] / R[numcells]);
	uc_left = U[0] - sqrt(GAMMA*P[0] / R[0]);
	uc_right = U[numcells - 5] - sqrt(GAMMA*P[numcells - 5] / R[numcells - 5]);
	printf("U-C left: %8.6lf\nU-C right: %8.6lf\nmiddle: %8.6lf\n", uc_left, uc_right, (uc_left + uc_right) / 2);
	system("pause");

#elif (PROBLEM == 2)
	double ro_right = 1.271413930046081;
	double ro_left = 1.0;
	double u_right = 0.292868067614595;
	double u_left = 0;

	D_analit = (ro_right * u_right - ro_left * u_left) / (ro_right - ro_left);
	//	printf("Analitical speed: %10.8lf\n\n", D_analit);
	//	system("pause");


#elif (PROBLEM == 8 || PROBLEM == 4)
	/*	printf("-----------OLD WAY----------\n");
	printf("\n--------first sw-------\n");
	printf("dens: %lf\n", gyugonio(st_P3, st_R3, st_P2));
	printf("sw speed: %lf\n", sw_speed(st_R3, st_R2, st_U3, st_U2));
	printf("u after sw: %lf\n", st_U2);
	printf("\n---------second sw---------\n");
	printf("dens: %lf\n", gyugonio(st_P2, st_R2, st_P1));
	printf("sw speed: %lf\n", sw_speed(st_R2, st_R1, st_U2, st_U1));


	printf("\n---------summa---------\n");
	printf("sw speed: %lf\n", sw_speed(st_R3, st_R1, st_U3, st_U1));
	D_analit = sw_speed(st_R3, st_R1, st_U3, st_U1);*/
	//D_analit = 3.48;

	printf("---------------NEW WAY-----------\n");
	printf("\n--------first sw-------\n");
	//printf("dens after sw: %lf\n", gyugonio(st_P3, st_R3, st_P2));
	double D_sw1 = sw_speed2(st_R3, st_U3, st_P3, gyugonio(st_P3, st_R3, st_P2), st_P2);
	printf("sw speed1 without u3: %lf\n", D_sw1);
	printf("rup after sw1: %lf, %lf, %lf\n", st_R2, st_U2, st_P2);
	system("pause");
	printf("\n--------second sw-------\n");
	//printf("dens after sw: %lf\n", st_R1);
	double D_sw2 = sw_speed2(st_R2, st_U2, st_P2, st_R1, st_P1);
	printf("sw speed2 without u3: %lf\n", D_sw2);
	printf("rup after sw2: %lf, %lf, %lf\n", st_R1, st_U1, st_P1);
	system("pause");
#endif
}

/**************************************************/

double* finite_difference(int numb, double *mas)
{
	int val = nmesh[numb];
	double *dif = new double[val - 2];
	int check = 0;
	long int p = 0;
	long int m = 0;
	long int zero = 0;


	for (int i = 0; i < val - 2; i++)
	{
		dif[i] = mas[i + 2] - 2 * mas[i + 1] + mas[i];
		if (i == 0)
		{
			if (dif[i] > 0) { check = 1; p++; continue; }
			if (dif[i] < 0) { check = 0; m++; continue; }
			if (dif[i] == 0) { check = 2; zero++; continue; }
		}

		if (dif[i] > 0)
		{
			if (check == 1) { p++; check = 1; }
			else {
				p++;
				if (check == 0) { if (m > 30) { printf("-: %d\n", m); } m = 0; }
				if (check == 2) { if (zero > 30) { printf("0: %d\n", zero); } zero = 0; }
				check = 1;
			}

		}
		if (dif[i] < 0)
		{
			if (check == 0) { m++; check = 0; }
			else {
				m++;
				if (check == 1) { if (p > 30) { printf("+: %d\n", p); } p = 0; }
				if (check == 2) { if (zero > 30) { printf("0: %d\n", zero); } zero = 0; }
				check = 0;
			}
		}

		if (dif[i] == 0)
		{
			if (check == 2) { zero++; check = 2; }
			else {
				zero++;
				if (check == 0) { if (m > 30) { printf("-: %d\n", m); } m = 0; }
				if (check == 1) { if (p > 30) { printf("+: %d\n", p); } p = 0; }
				check = 2;
			}
		}

		if (i == val - 3)
		{
			if (check == 0) printf("-: %d\n", m);
			if (check == 1) printf("+: %d\n", p);
			if (check == 2) printf("0: %d\n", zero);
		}

	}

	return dif; // возвращаем указатель на массив!!! его имя!
}

void file_exact_diff(int numcells, double *exact_R, double *exact_U, double *exact_P, double *exact_RE, double *exact_S, double *diff_R, double *diff_U, double *diff_P, double time)
{
}

void runge(double *massiv, int lda, int numb)
{
	double value[NUM_ITER - 2]; // without boundary points
	for (int i = 0; i < NUM_ITER - 2; i++)
	{
		//value[i] = log(fabs((fabs(massiv[i]) - fabs(massiv[i + 1])) / (fabs(massiv[i + 1]) - fabs(massiv[i + 2])))) / log(double(GRID));
		value[i] = log(fabs((massiv[numb + lda * i] - massiv[numb + lda * (i + 1)]) / (massiv[numb + lda * (i + 1)] - massiv[numb + lda * (i + 2)]))) / log(double(GRID));

		printf("Precision order for iteration double %d: %lf\n", i + 1, value[i]);
	}
}

void analitical_writing_into_file(int numcells, double* R_D, double*R_U, double*R_P, double timer)
{
	FILE* fout;
	double x, x_NC;
	double dx = LENGTH / double(numcells);
	char name1[255];

	double D_analit;
	double ro_right, ro_left, u_right, u_left;

	ro_right = 1.271413930046081;
	ro_left = 1.0;
	u_right = 0.292868067614595;
	u_left = 0;

	D_analit = (ro_right * u_right - ro_left * u_left) / (ro_right - ro_left);

#ifndef NC
	sprintf(name1, "workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_analit_%6.4lf.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, timer);
#else
	sprintf(name1, "workspace/%03d/NC_N%03d_P%1d_SLV%1d_TERM%.0lf_analit_%6.4lf.dat", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, timer);
#endif

	fout = fopen(name1, "w");

	for (int i = 0; i < numcells; i++)
	{

#ifndef NC
		x = i*dx + 0.5*dx;
#else
#ifdef NC2
		x_NC = (i*dx + 0.5*dx - DISC_POINT - D_analit*timer) / dx;
#else
		x_NC = i*dx + 0.5*dx - DISC_POINT - D_analit*timer;
#endif

#endif


		/**********************************
		| 1 |    2    |   3   |     4    |
		| x | density | speed | pressure |
		***********************************/
#ifndef NC
		fprintf(fout, "%9.6lf %lf %lf %lf \n", x, R_D[i], R_U[i], R_P[i]);
#else
		fprintf(fout, "%9.6lf %lf %lf %lf \n", x_NC, R_D[i], R_U[i], R_P[i]);
#endif
	}
	fclose(fout);
}

