#include "definitions.h"
#include "support.h"

/****************************************************
Source file contains functions responsible
for printing the results of output
using GNUPLOT application.

To enable automatic output, please set #define PRINT in the definitions.h

The output directory must consists of special folders.
The structure is the following:

If you use Windows OS with Microsoft Visual Studio IDE, build the folder 'workspace' in the folder with .vcxproj file.
If you use Linux OS, build the folder 'workspace' in the folder with a.out file.

The 'workspace' folder must have the structure:
1) The folders, corresponding to the used grid points: 100, 300, 900, ..., 24300,...
2) Each of this folders should have five subfolders: D, P, S, R, U - for the gas values.
3) Each of the previous folder should contain the subfolder with the name of the problem: P_0, P_1, P_2, ..., P_20, ...

****************************************/

void gnuplot_n_smooth(int numb)
{
	FILE* plot;
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / double(N_smooth);
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
		time_control[i] = (i + 1)*k_step;

	for (int i = 0; i < 3; i++)
	{
		plot = fopen("N_smooth.plt", "w");
		fprintf(plot, "reset\nclear\n\n###################################\nset term png font \"Times-Roman, 16\"\n\n##################################\n\n");
		fprintf(plot, "set xrange[0:%f]\n\n", LENGTH); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
		if (PROBLEM == 0) fprintf(plot, "set yrange[%4.3f:%4.3f]\n\n", left_SW[i], right_SW[i]);
		if (PROBLEM == 1) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW[i], right_RW[i]);
		if (PROBLEM == 2) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_ST[i], right_ST[i]);
		if (PROBLEM == 3) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_SW_cs[i], right_SW_cs[i]);
		if (PROBLEM == 4) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_SW[i], right_SW[i]);
		if (PROBLEM == 5) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_2RR[i], right_2RR[i]);
		if (PROBLEM == 6) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW_SW[i], right_RW_SW[i]);
		if (PROBLEM == 7 || PROBLEM == 8) fprintf(plot, "set yrange[%4.3f:%4.3f]\n\n", left_RW_RW[i], right_RW_RW[i]);
		//if (PROBLEM == 8) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW1_RW1[i], right_RW1_RW1[i]);
		if (PROBLEM == 9) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_sodd[i], right_sodd[i]);
		if (i == 0) fprintf(plot, "set ylabel \"density\"\n");
		if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\n");
		if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\n");
		fprintf(plot, "set xlabel \"x\"\n");


		for (int k = 0; k < N_smooth; k++)
		{
#ifdef TIME
			if ((time_control[k] - 0.0015) < k_step || (time_control[k] - 0.006) < k_step || (time_control[k] - 0.025) < k_step || (time_control[k] - 0.1) < k_step)
			{
#endif
				fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/LP%03d_P%1d_%6.4lf_%c.png'\n", nmesh[numb], prop[i], PROBLEM, nmesh[numb], PROBLEM, (k + 1)*k_step, prop[i]);
				//fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/SZ%03d_P%1d_%6.4lf_%c.png'\n", nmesh[numb], prop[i], PROBLEM, nmesh[numb], PROBLEM, (k + 1)*k_step, prop[i]);
				//fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4f.dat' using 1:%d w l lw 2 notitle\n\n", nmesh[numb], nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (k + 1)* k_step, i + 2);
				fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4f.dat' using 1:%d w linespoints pt 7 notitle\n\n", nmesh[numb], nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (k + 1)* k_step, i + 2);
#ifdef TIME
			}
#endif
		}
		fprintf(plot, "\nquit\n");
		fclose(plot);
		system("N_smooth.plt");
	}

	return;
}

void gnuplot_n_smooth2(int numcells, int sw1[3][N_smooth], int sw2[3][N_smooth], int sw3[3][N_smooth])
{
	FILE* plot;
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
		time_control[i] = (i + 1)*k_step;

	//RUP
	double left[5];
	double right[5];

	int small = 0;

	left[0] = st_th_R2;
	left[1] = st_th_U2;
	left[2] = st_th_P2;
	left[3] = 0;
	left[4] = sqrt(GAMMA*st_th_P2 / st_th_R2) - 1;

	right[0] = st_th_R1;
	right[1] = st_th_U1;
	right[2] = st_th_P1;
	right[3] = 0;
	right[4] = sqrt(GAMMA*st_th_P1 / st_th_R1) + 1;

	for (int i = 0; i < 6; i++)
	{
		if (i == 5)
		{
			// SMALL OR BIG
			small = 0;
		}
		if (i == 3) continue;  // entropy - 4 , diff in entropy - 5
		plot = fopen("N_smooth2.plt", "w");
		for (int k = 0; k < N_smooth; k++)
		{
			
			fprintf(plot, "\nreset\nclear\n\n###################################\nset term png font \"Times-Roman, 16\"\n\n##################################\n\n");

			fprintf(plot, "set xrange[0.0:%f]\n\n", LENGTH); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
			fprintf(plot, "set xlabel \"x\"\n");

			// Boundary conditions
			set_bound(plot, i);

			if (i == 0) fprintf(plot, "set ylabel \"density\"\n");
			if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\n");
			if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\n");
			if (i == 4)	fprintf(plot, "set ylabel \"entropy\"\n");
			if (i == 5)	fprintf(plot, "set ylabel \"entropy difference\"\n");
#ifdef SW_POINTS_PRINT
			if (i == 0) fprintf(plot, "f(x) = %lf\ng(x) = %lf\n", initial_density(0.03) - DELTA, initial_density(0.95) + DELTA);
			if (i == 1)	fprintf(plot, "f(x) = %lf\ng(x) = %lf\n", initial_velocity(0.03) - DELTA, initial_velocity(0.95) + DELTA);
			if (i == 2)	fprintf(plot, "f(x) = %lf\ng(x) = %lf\n", initial_pressure(0.03) - DELTA, initial_pressure(0.95) + DELTA);

#if(PROBLEM == 8 || PROBLEM == 4)
			if (time_control[k] < 0.26)
			{
				if (i == 0) fprintf(plot, "\nset label \"N_{sw1} = %d\" at 0.1, 0.2\n\n", sw1[0][k]);
				if (i == 1) fprintf(plot, "\nset label \"N_{sw1} = %d\" at 0.1, 0.2\n\n", sw1[1][k]);
				if (i == 2) fprintf(plot, "\nset label \"N_{sw1} = %d\" at 0.1, 0.2\n\n", sw1[2][k]);

				if (i == 0) fprintf(plot, "\nset label \"N_{sw2} = %d\" at 0.1, 0.0\n\n", sw2[0][k]);
				if (i == 1) fprintf(plot, "\nset label \"N_{sw2} = %d\" at 0.1, 0.0\n\n", sw2[1][k]);
				if (i == 2) fprintf(plot, "\nset label \"N_{sw2} = %d\" at 0.1, 0.0\n\n", sw2[2][k]);
			}
			else
			{
				if (i == 0) fprintf(plot, "\nset label \"N_{sw3} = %d\" at 0.1, -0.2\n\n", sw3[0][k]);
				if (i == 1) fprintf(plot, "\nset label \"N_{sw3} = %d\" at 0.1, -0.2\n\n", sw3[1][k]);
				if (i == 2) fprintf(plot, "\nset label \"N_{sw3} = %d\" at 0.1, -0.2\n\n", sw3[2][k]);
			}
#endif
#endif		

#ifdef NC2
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/%c_NC2_N%03d_P%1d_%6.4lf_%c_NC.png'\n", numcells, prop[i], PROBLEM, (char)TYPE, numcells, PROBLEM, (k + 1)*k_step, prop[i]);
#else
#ifndef P4_ONE_WAVE
			if (small == 1)
			{
				fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/small_%c_W%03d_P%1d_%6.4lf_%c.png'\n", numcells, prop[i], PROBLEM, (char)TYPE, numcells, PROBLEM, (k + 1)*k_step, prop[i]);
			}
			else {
				fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/%c_W%03d_P%1d_%6.4lf_%c.png'\n", numcells, prop[i], PROBLEM, (char)TYPE, numcells, PROBLEM, (k + 1)*k_step, prop[i]);
			}
#else
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/%c_%03d_P%1d_%6.4lf_%c.png'\n", numcells, prop[i], PROBLEM, (char)TYPE, numcells, PROBLEM, (k + 1)*k_step, prop[i]);
#endif
#endif

			fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4f.dat' using 1:%d w linespoints pt 7 lt rgb 'purple' notitle\n\n", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (char)TYPE, (k + 1)* k_step, i + 2);

		}
		small = 0;
		fclose(plot);
		system("N_smooth2.plt");
	}


	return;
}

void gnuplot_all_iter_one_time(int run[NUM_ITER], int numb, double time)
{
	FILE* plot;
	char FileName[255];
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	for (int i = 0; i < N_smooth; i++)
	{
		time_control[i] = (i + 1)*k_step;
		if (fabs(time_control[i] - time) < EPS) time = time_control[i]; // do not forget to set up an absolute value
	}
	for (int i = 0; i < 6; i++) // for each parameter of gas
	{
		if (i == 3) continue;
		plot = fopen("all_iter.plt", "w");
		set_bound(plot, i);
		fprintf(plot, "set xlabel \"x\"\n");
		if (i == 0) fprintf(plot, "set ylabel \"density\"\n");
		if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\n");
		if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\n");
		if (i == 4)	fprintf(plot, "set ylabel \"entropy\"\n");
		if (i == 5)	fprintf(plot, "set ylabel \"entropy difference\"\n");
		fprintf(plot, "set output 'workspace/all_iter_one_time%6.4f_P%1d_%c_%c.png' \n\n", time, PROBLEM, (char)TYPE, prop[i]);
		fprintf(plot, "plot ");

		for (int j = 0; j <= numb; j++)
		{
			if (run[j] != 0)
			{
				int numcells = 100 * pow(3, j);
				fprintf(plot, "'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4f.dat' using 1:%d w linespoints pt 7 title \"N = %03d\"", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (char)TYPE, time, i + 2, numcells);
				if (j < numb) fprintf(plot, ", ");
			}

		}
		fclose(plot);
		system("all_iter.plt");
	}
}

void set_bound(FILE *plot, int i)
{
		fprintf(plot, "set term png font \"Times - Roman, 16\"\n\n");
		if (PROBLEM == 0 || PROBLEM == 18) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_SW[i], right_SW[i]);
		if (PROBLEM == 1) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_RW[i], right_RW[i]);
		if (PROBLEM == 2 || PROBLEM == 17) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_ST[i], right_ST[i]);
		if (PROBLEM == 3)  fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", -1.0, 1.0);
		if (PROBLEM == 8 || PROBLEM == 20)  fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_RW_RW[i], right_RW_RW[i]);
		if (PROBLEM == 5) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_2RR[i], right_2RR[i]);
		if (PROBLEM == 19) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_19[i], right_19[i]);
		if (PROBLEM == 7) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_7[i], right_7[i]);
		if (PROBLEM == 4) fprintf(plot, "set yrange[%9.8f:%9.8f]\n\n", left_SW_SW[i], right_SW_SW[i]);
	
}

void gnuplot_n_smooth3(int numcells)
{
	FILE* plot;
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
		time_control[i] = (i + 1)*k_step;

	//RUP
	double left[5];
	double right[5];

	left[0] = st_th_R2;
	left[1] = st_th_U2;
	left[2] = st_th_P2;
	left[3] = 0;
	//left[4] = sqrt(GAMMA*st_th_P2 / st_th_R2) - 1;
	left[4] = 0.8;

	right[0] = st_th_R1;
	right[1] = st_th_U1;
	right[2] = st_th_P1;
	right[3] = 0;
	//right[4] = sqrt(GAMMA*st_th_P1 / st_th_R1) + 1;
	right[4] = 1.3;

	for (int i = 0; i < 5; i++)
	{
		if (i == 3) continue;  // entropy - 4
		plot = fopen("N_smooth3.plt", "w");
		for (int k = 0; k < N_smooth; k++)
		{

			fprintf(plot, "\nreset\nclear\n\n###################################\nset term png font \"Times-Roman, 16\"\n\n##################################\n\n");

			fprintf(plot, "set xrange[0.0:1.0]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
			fprintf(plot, "set xlabel \"x\"\n");

			if (PROBLEM == 0) fprintf(plot, "set yrange[%4.3f:%4.3f]\n\n", left_SW[i], right_SW[i]);
			if (PROBLEM == 1) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW[i], right_RW[i]);
			if (PROBLEM == 8)  fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW_RW[i], right_RW_RW[i]);
			if (PROBLEM == 10) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left[i], right[i]);

#ifdef SW_POINTS_PRINT
			if (i == 0) fprintf(plot, "set ylabel \"density\"\nf(x) = %lf\ng(x) = %lf\n", initial_density(0.03) - DELTA, initial_density(0.95) + DELTA);
			if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\nf(x) = %lf\ng(x) = %lf\n", initial_velocity(0.03) - DELTA, initial_velocity(0.95) + DELTA);
			if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\nf(x) = %lf\ng(x) = %lf\n", initial_pressure(0.03) - DELTA, initial_pressure(0.95) + DELTA);
			if (i == 4)	fprintf(plot, "set ylabel \"entropy\"\n");

#endif
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/W%03d_P%1d_%6.4lf_%c.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, (k + 1)*k_step, prop[i]);

			fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4f.dat' using 1:%d w linespoints pt 7 title 'exact RP', \
						  'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4f.dat' using 1:%d w linespoints pt 7 title 'linear RP' \n\n",
				numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, 'E', (k + 1)* k_step, i + 2, numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, 'L', (k + 1)* k_step, i + 2);

		}
		fclose(plot);
		system("N_smooth3.plt");
	}


	return;
}

void gnuplot_one_iteration(int numcells) //together with analitical solution
{
	FILE *plot;
	double iter;

	iter = log(float(numcells / 100)) / log(3.0) + 1;

	for (int i = 0; i < 3; i++) //итерации по характеристикам газа
	{
		plot = fopen("Plot_new.plt", "w");
		fprintf(plot, "reset\n\
					   clear\n\n\
					   ###################################\n\
						set term png font \"Times-Roman, 16\"\n\n\
						##################################\n\n");
#ifdef NC
#if (PROBLEM==0)
		fprintf(plot, "set xrange[-0.3:0.7]\n\n");
		fprintf(plot, "set yrange[%3.2f:%9.7f]\n\
					   set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_NORP_%c_NC.png'\n\
					   plot 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1 : %d w l lw 3 notitle\n\n", left_SW[i], right_SW[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, prop[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
#elif (PROBLEM==1)
		fprintf(plot, "set xrange[0:5]\n\n");
		fprintf(plot, "set yrange[%3.2f:%9.7f]\n\
					   set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_NORP_%c_NC.png'\n\
					   plot 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1 : %d w l lw 3 notitle", left_RW[i], right_RW[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, prop[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
#elif(PROBLEM==2)
		fprintf(plot, "set xrange[0:5]\n\n");
		fprintf(plot, "set yrange[%3.2f:%9.7f]\n\
					   set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_NORP_%c_NC.png'\n\
					   plot 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1 : %d w l lw 3 notitle", left_RW[i], right_RW[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, prop[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);

#endif
#else
		fprintf(plot, "set xrange[0:1]\n\n");
#if (PROBLEM==0)
		fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_SW[i], right_SW[i]);
#elif (PROBLEM==1)
		fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_RW[i], right_RW[i]);
		analitical_RW(plot, initial_pressure(0.05), initial_density(0.05), initial_velocity(0.05), initial_pressure(0.2), initial_density(0.2), initial_velocity(0.2), iter);  // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
#elif (PROBLEM==2)
		fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_RP[i], right_RP[i]);
#endif

		fprintf(plot, "set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_NORP_%c.png'\nplot 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1 : %d w l lw 3 notitle", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, prop[i], numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);

#if (PROBLEM==0) //analitical solution
		switch (i)
		{
		case 0:
			fprintf(plot, ", 'data_dl.txt' w l lw 3 lt 2 notitle, 'data_dm.txt' w l lw 3 lt 2 notitle, 'data_dr.txt' w l lw 3 lt 2 notitle\n\n"); break; // если у нас волна разрежения, то вместе с приближенной волной рисуем АНАЛИТИЧЕСКУЮ волну с ОСТРЫМИ концами
		case 1:
			fprintf(plot, ", 'data_ul.txt' w l lw 3 lt 2 notitle, 'data_um.txt' w l lw 3 lt 2 notitle, 'data_ur.txt' w l lw 3 lt 2 notitle\n\n"); break;
		case 2:
			fprintf(plot, ", 'data_pl.txt' w l lw 3 lt 2 notitle, 'data_pm.txt' w l lw 3 lt 2 notitle, 'data_pr.txt' w l lw 3 lt 2 notitle\n\n"); break;
		default:
			break;
		}

#elif (PROBLEM==1)
		switch (i)
		{
		case 0: fprintf(plot, ", ro(x) w l lw 3 lt 2 notitle, 'data_dl.txt' w l lw 3 lt 2 notitle, 'data_dr.txt' w l lw 3 lt 2 notitle\n\n"); break; // если у нас волна разрежения, то вместе с приближенной волной рисуем АНАЛИТИЧЕСКУЮ волну с ОСТРЫМИ концами
		case 1: fprintf(plot, ", u(x) w l lw 3 lt 2 notitle, 'data_ul.txt' w l lw 3 lt 2 notitle, 'data_ur.txt' w l lw 3 lt 2 notitle\n\n"); break;
		case 2: fprintf(plot, ", p(x) w l lw 3 lt 2 notitle, 'data_pl.txt' w l lw 3 lt 2 notitle, 'data_pr.txt' w l lw 3 lt 2 notitle\n\n"); break;
		default: break;
		}

#endif
#endif

		fclose(plot);
		system("Plot_new.plt");
	}

}

void gnuplot_RW_DIFF(int numcells)
{
	FILE *plot;
	double iter;

	iter = log(float(numcells / 100)) / log(3.0) + 1;

	plot = fopen("Plot_RW_DIFF.plt", "w");
	fprintf(plot, "reset\n\
				   clear\n\n\
				   ##################################\n\
  				   set term png font \"Times-Roman, 16\"\n\n\
				   ##################################\n\n");
#ifndef NC
	fprintf(plot, "set xrange[0:1]\n\n");
#else
	//fprintf(plot, "set xrange[0.74:1.26]\n\n");
	fprintf(plot, "set xrange[0:2]\n\n");
#endif
	fprintf(plot, "set yrange[-0.2:0.2]\n\
				   set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_RW_DIFF_U.png'\n\
												  plot 'N%04d_RW_difference.dat' using 1 : 2 w l lw 3 notitle\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, numcells);
	fprintf(plot, "set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_RW_DIFF_P.png'\n\
				  				  				  plot 'N%04d_RW_difference.dat' using 1 : 3 w l lw 3 notitle\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, numcells);
	fprintf(plot, "set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_RW_DIFF_R.png'\n\
				  				  				  plot 'N%04d_RW_difference.dat' using 1 : 4 w l lw 3 notitle\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, numcells);

	fclose(plot);
	system("Plot_RW_DIFF.plt");
}

void gnuplot_RW_NUM_ANALITIC(int numcells)
{
	FILE *plot;

	plot = fopen("Plot_RW_NUM_ANALITIC.plt", "w");
	fprintf(plot, "reset\n\
				  clear\n\n\
				 ###################################\n\
				 set term png font \"Times-Roman, 16\"\n\n\
				 ##################################\n\n");
#ifndef NC
	fprintf(plot, "set xrange[0:1]\n\n");
#else
	fprintf(plot, "set xrange[0.9:1.1]\n\n");
#endif
	fprintf(plot, "set yrange[-0.35:0.1]\n\
				  				  				  set output 'N%04d_P%1d_NUM_ANALITIC_U.png'\n\
												  								  				  plot 'N%04d_RW_NUM_ANALITIC.dat' using 1 : 2 w l lw 2 notitle, 'N%04d_RW_NUM_ANALITIC.dat' using 1 : 3 w l lw 2 notitle\n\n", numcells, PROBLEM, numcells, numcells);
	fprintf(plot, "set yrange[1.3:2.1]\nset output 'N%04d_P%1d_NUM_ANALITIC_P.png'\n\
				  				  				  plot 'N%04d_RW_NUM_ANALITIC.dat' using 1 : 4 w l lw 2 notitle, 'N%04d_RW_NUM_ANALITIC.dat' using 1 : 5 w l lw 2 notitle\n\n", numcells, PROBLEM, numcells, numcells);
	fprintf(plot, "set yrange[1.5:2.1]\nset output 'N%04d_P%1d_NUM_ANALITIC_R.png'\n\
				  				  				  plot 'N%04d_RW_NUM_ANALITIC.dat' using 1 : 6 w l lw 2 notitle, 'N%04d_RW_NUM_ANALITIC.dat' using 1 : 7 w l lw 2 notitle\n\n", numcells, PROBLEM, numcells, numcells);

	fclose(plot);
	system("Plot_RW_NUM_ANALITIC.plt");
}

void gnuplot_P_PLUS_PG(int numcells)
{

	FILE *plot;

	plot = fopen("Plot_P_PLUS_PG.plt", "w");
	fprintf(plot, "reset\n\
				   clear\n\n\
				   ###################################\n\
				   set term png font \"Times-Roman, 16\"\n\n\
				   #################################\n\n");
	fprintf(plot, "set xrange[-0.1:0.3]\n\n");
	fprintf(plot, "set yrange[-0.1:1.5]\n\
				   set output 'N%04d_P%1d_SLV%1d_TERM%.0lf_P_PLUS_PG.png'\n\
				   plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_P_PLUS_PG_NC.dat' using 1 : 2 w l lw 3 notitle, 'N%04d_P%1d_SLV%1d_TERM%.0lf_P_PLUS_PG_NC.dat' using 1 : 3 w l lw 3 notitle\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);

	fclose(plot);
	system("Plot_P_PLUS_PG.plt");

}

void gnuplot_all_iterations_NC(int numb)
{
	FILE *plot;
	plot = fopen("Plot_new_all_iterations.plt", "w");
	fprintf(plot, "reset\n\
				   clear\n\n\
				   set term png font \"Times-Roman, 16\"\n\n\
				   ##################################\n\n");
	for (int i = 0; i < 3; i++)  // choosing gas property
	{
#if (PROBLEM==0)
		fprintf(plot, "set xrange[-0.3:0.7]\n");
		fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_SW[i], right_SW[i]);
#elif (PROBLEM==1)
#ifdef RW_NUM_ANALITICAL
		fprintf(plot, "set xrange[0:2]\n\n");
		fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_RW[i], right_RW[i]);
#else
		fprintf(plot, "set xrange[0:5]\n\n");
		fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_RW[i], right_RW[i]);
#endif

#endif

		fprintf(plot, "set output 'NSUM_P%1d_SLV%1d_TERM%.0lf_NORP_%c_NC.png'\n", PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, prop[i]);

		fprintf(plot, "plot ");
		for (int j = 0; j < numb; j++)
		{
			fprintf(plot, "'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1 : %d w l lw 2 notitle, ", nmesh[j], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
			if (j == numb - 1) fprintf(plot, "'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1 : %d w l lw 2 notitle\n\n", nmesh[j], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
		}

	}
	fclose(plot);
	system("Plot_new_all_iterations.plt");
}

void gnuplot_one_it_NC()
{
	float const bound = 0.1;
	FILE *plot;
	plot = fopen("Plot_one_it_NC.plt", "w");
	fprintf(plot, "reset\n\
				  				  				  				  clear\n\n\
																  												  								  				  set term png font \"Times-Roman, 16\"\n\n\
																																								  																								  												  					 ##################################\n\n");
	double step = time_max_array[PROBLEM] / N_smooth;
	for (int k = 0; k < NUM_ITER; k++)
	{
		for (int i = 0; i < 3; i++)  // choosing gas property
		{
#if (PROBLEM==0)
			fprintf(plot, "set xrange[-0.1:0.1]\n");
			fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_SW[i], right_SW[i]);
#elif (PROBLEM==1)
#ifdef RW_NUM_ANALITICAL
			fprintf(plot, "set xrange[0:0.3]\n\n");
			fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_RW[i], right_RW[i]);
#else
			fprintf(plot, "set xrange[0:0.3]\n\n");
			fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_RW[i], right_RW[i]);
#endif
#elif PROBLEM==7
			fprintf(plot, "set xrange[-10:10]\n");
			fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW_RW[i], right_RW_RW[i]);
#endif
			fprintf(plot, "set xlabel \"X\"\n");
			if (i == 0) fprintf(plot, "set ylabel \"density\"\n");
			if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\n");
			if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\n");

			//		if (k == 0)
			//		{
			//			fprintf(plot, "set xrange[-%f:%f]\n", bound / pow(3.0, k), bound / pow(3.0, k));
#ifdef alpha
			fprintf(plot, "set output 'workspace/%03d/%c/NSUM%03d_P%1d_%c_one_it_NC_C%3.1f_alpha%3.2f.png'\n", nmesh[k], prop[i], nmesh[k], PROBLEM, prop[i], C1, alpha);
#else
			fprintf(plot, "set output 'workspace/%03d/%c/NC_NSUM%03d_P%1d_%c_one_it.png'\n", nmesh[k], prop[i], nmesh[k], PROBLEM, prop[i]);
#endif
			fprintf(plot, "plot ");


			for (int j = 0; j < N_smooth; j++)
			{
#if (PROBLEM == 7)
				if ((j + 1)*step >= 4.0) fprintf(plot, "'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4lf_NC.dat' using 1 : %d w l lw 2 notitle, ", nmesh[k], nmesh[k], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (j + 1)*step, i + 2);
#else
				fprintf(plot, "'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4lf_NC.dat' using 1 : %d w l lw 2 notitle, ", nmesh[k], nmesh[k], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (j + 1)*step, i + 2);
#endif
				if (j == N_smooth - 1) fprintf(plot, "'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4lf_NC.dat' using 1 : %d w l lw 2 notitle\n\n", nmesh[k], nmesh[k], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (j + 1)*step, i + 2);
			}
			//	}



		}
	}
	fclose(plot);
	system("Plot_one_it_NC.plt");
	return;
}

void gnuplot_conservative(int numb)
{
	FILE *plot;

	for (int i = 0; i < 3; i++) //итерации по характеристикам газа
	{
		plot = fopen("Plot_conserv.plt", "w");
		fprintf(plot, "reset\n\
					  					  clear\n\n\
										  					   ###################################\n\
															   					   set term png font \"Times-Roman, 16\"\n\n\
																				   					   ##################################\n\n");

		fprintf(plot, "set xrange[0:1]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения

		for (int j = 0; j < 4; j++)
		{

			fprintf(plot, "set output 'N%04d_P%1d_T%c_5_%d.png'\n", nmesh[numb], PROBLEM, dip[i], j + 1);
			fprintf(plot, "plot 'N%04d_P%1d_T%c.dat' using 1:%d w l notitle\n\n", nmesh[numb], PROBLEM, dip[i], j + 2);
		}


		fclose(plot);
		system("Plot_conserv.plt");
	}
}

void gnuplot_five_t_steps(int numb)
{
	FILE *plot;

	for (int i = 0; i < 3; i++) //итерации по характеристикам газа
	{
		plot = fopen("Plot_five_it.plt", "w");
		fprintf(plot, "reset\n\
					  					  					  clear\n\n\
															  										  						###################################\n\
																																																					set term png font \"Times-Roman, 16\"\n\n\
																																																																																	##################################\n\n");

		fprintf(plot, "set xrange[0:1]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
		fprintf(plot, "set yrange[-1:1]\n\n");

		for (int j = 0; j < N_smooth; j++) // шги по времени
		{
			fprintf(plot, "set output 'N%04d_P%1d_%4.2lf_%c.png'\n", nmesh[numb], PROBLEM, (j + 1)*0.05, prop[i]);
			fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_%4.2f.dat' using 1:%d w l notitle\n\n", nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (j + 1)* 0.05, i + 2);

		}


		fclose(plot);
		system("Plot_five_it.plt");
	}

}

void gnuplot_n_smooth_NC(int numb)
{
	FILE* plot;
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
		time_control[i] = (i + 1)*k_step;

#ifdef TIME
	if (numb == TIME)
	{
#endif
		for (int i = 0; i < 3; i++)
		{
			plot = fopen("N_smooth_NC.plt", "w");
			fprintf(plot, "reset\nclear\n\n###################################\nset term png font \"Times-Roman, 16\"\n\n##################################\n\n");
			//	fprintf(plot, "set xrange[-0.5:0.5]\n\n");
#ifdef NC2
			fprintf(plot, "set xrange[-15:15]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
			fprintf(plot, "set xlabel \"(x-Dt-x0)/h\"\n");
#endif
			if (PROBLEM == 0) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_SW[i], right_SW[i]);
			if (PROBLEM == 1) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW[i], right_RW[i]);

#ifdef SW_POINTS_PRINT
			if (i == 0) fprintf(plot, "set ylabel \"density\"\nf(x) = %lf\ng(x) = %lf\n", initial_density(0.03) - DELTA, initial_density(0.95) + DELTA);
			if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\nf(x) = %lf\ng(x) = %lf\n", initial_velocity(0.03) - DELTA, initial_velocity(0.95) + DELTA);
			if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\nf(x) = %lf\ng(x) = %lf\n", initial_pressure(0.03) - DELTA, initial_pressure(0.95) + DELTA);
#endif
			for (int k = 0; k < N_smooth; k++)
			{
#ifdef TIME
				if (fabs(time_control[k] - 0.0015) < k_step || fabs(time_control[k] - 0.0035) < k_step || fabs(time_control[k] - 0.015) < k_step || fabs(time_control[k] - 0.08) < k_step)
				{
#endif
#ifdef NC2
#ifdef SW_POINTS_PRINT
					fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/pr_swp_N%03d_P%1d_%6.4lf_%c_NC.png'\n", nmesh[numb], prop[i], PROBLEM, nmesh[numb], PROBLEM, (k + 1)*k_step, prop[i]);
#else
					fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/swp_N%03d_P%1d_%6.4lf_%c_NC.png'\n", nmesh[numb], prop[i], PROBLEM, nmesh[numb], PROBLEM, (k + 1)*k_step, prop[i]);
#endif
					fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4f_NC.dat' using 1:%d w linespoints pt 7 notitle, f(x) lt rgb 'green' notitle, g(x) lt rgb 'green' notitle\n\n", nmesh[numb], nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (k + 1)* k_step, i + 2);
#else
					fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/N%03d_P%1d_%6.4lf_%c_NC.png'\n", nmesh[numb], prop[i], PROBLEM, nmesh[numb], PROBLEM, (k + 1)*k_step, prop[i]);
					fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4f_NC.dat' using 1:%d w l notitle\n\n", nmesh[numb], nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (k + 1)* k_step, i + 2);
#endif
					//	fprintf(plot, "set output 'N%04d_P%1d_%6.4lf_%c_NC.png'\n", nmesh[numb], PROBLEM, time_control[k], prop[i]);
					//	fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_%6.4f_NC.dat' using 1:%d w linespoints pt 7 notitle, f(x) lt rgb 'green' notitle, g(x) lt rgb 'green' notitle\n\n", nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 2);
#ifdef TIME
				}
#endif
			}

			fclose(plot);
			system("N_smooth_NC.plt");
		}
#ifdef TIME
	}
#endif

	return;
}

void gnuplot_n_smooth_NC2(int numcells, int* n_r, int* n_u, int* n_p)
{
	FILE* plot;
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
		time_control[i] = (i + 1)*k_step;


	for (int i = 0; i < 3; i++)
	{
		plot = fopen("N_smooth_NC.plt", "w");
		for (int k = 0; k < N_smooth; k++)
		{

			fprintf(plot, "\nreset\nclear\n\n###################################\nset term png font \"Times-Roman, 16\"\n\n##################################\n\n");

			fprintf(plot, "set xrange[-15:15]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
			fprintf(plot, "set xlabel \"(x-Dt-x0)/h\"\n");

			if (PROBLEM == 0) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_SW[i], right_SW[i]);
			if (PROBLEM == 1) fprintf(plot, "set yrange[%3.2f:%3.2f]\n\n", left_RW[i], right_RW[i]);

#ifdef SW_POINTS_PRINT
			if (i == 0) fprintf(plot, "set ylabel \"density\"\nf(x) = %lf\ng(x) = %lf\n", initial_density(0.05) - DELTA, initial_density(0.95) + DELTA);
			if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\nf(x) = %lf\ng(x) = %lf\n", initial_velocity(0.05) - DELTA, initial_velocity(0.95) + DELTA);
			if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\nf(x) = %lf\ng(x) = %lf\n", initial_pressure(0.05) - DELTA, initial_pressure(0.95) + DELTA);

#endif
#ifdef SW_POINTS_PRINT
			if (i == 0) fprintf(plot, "\nset label \"N_{sw} = %d\" at 7, 1.2\n\n", n_r[k]);
			if (i == 1) fprintf(plot, "\nset label \"N_{sw} = %d\" at 7, 0.25\n\n", n_u[k]);
			if (i == 2) fprintf(plot, "\nset label \"N_{sw} = %d\" at 7, 1.3\n\n", n_p[k]);

			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/pr_NC2_N%03d_P%1d_%6.4lf_%c_NC.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, (k + 1)*k_step, prop[i]);
#else
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/NC2_N%03d_P%1d_%6.4lf_%c_NC.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, (k + 1)*k_step, prop[i]);
#endif
			fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4f_NC.dat' using 1:%d w linespoints pt 7 notitle, f(x) lt rgb 'green' notitle, g(x) lt rgb 'green' notitle\n\n", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, (k + 1)* k_step, i + 2);

			//	fprintf(plot, "set output 'N%04d_P%1d_%6.4lf_%c_NC.png'\n", nmesh[numb], PROBLEM, time_control[k], prop[i]);
			//	fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_%6.4f_NC.dat' using 1:%d w linespoints pt 7 notitle, f(x) lt rgb 'green' notitle, g(x) lt rgb 'green' notitle\n\n", nmesh[numb], PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 2);
		}


		fclose(plot);
		system("N_smooth_NC.plt");
	}


	return;
}

void gnuplot_analitical_riemann(int numcells, double* R, double*U, double*P, double* R_D, double*R_U, double*R_P)
{
	FILE* fout;
	FILE* plot;
	char name[255];
	double dx = LENGTH / double(numcells);
	double x;

	sprintf(name, "N%04d_P%1d_SLV%1d_TERM%.0lf_riemann_analit.dat", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);
	fout = fopen(name, "w");


	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;

		/**********************************
		| 1 |    2    |   3   |     4    |
		| x | density | speed | pressure |
		***********************************/

		fprintf(fout, "%9.6lf %lf %lf %lf \n", x, R_D[i], R_U[i], R_P[i]);
	}

	/*----------------------------------------*/
	sprintf(name, "N%04d_P%1d_riemann_analit.plt", numcells, PROBLEM);
	plot = fopen(name, "w");

	fprintf(plot, "reset\n\
				  				  clear\n\n\
								  				   ###################################\n\
												   				   set term png font \"Times-Roman, 16\"\n\n\
																   				    ##################################\n\n");
	fprintf(plot, "set xlabel \"X\"\n");
	fprintf(plot, "set xrange[0:1]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения

	for (int i = 0; i < 3; i++) //итерации по характеристикам газа
	{
		if (i == 1) fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_SW[i], right_SW[i]);
		else fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_SW[i], right_RW[i]);
		fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/N%04d_P%1d_%c_analit_riemann.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, dip[i]);
		if (i == 0) fprintf(plot, "set ylabel \"density\"\n");
		if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\n");
		if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\n");
#ifdef RP
		fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_riemann_analit.dat' using 1:%d w l lw 2 title \"exact\", 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w linespoints pt 4 lt rgb 'blue' title \"num\", 'RP_N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w linespoints pt 7 lt rgb 'green' title \"num_rp\" \n\n",
			numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
#else
		fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_riemann_analit.dat' using 1:%d w l lw 2 title \"exact\", 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w linespoints pt 7 title \"numerical\"\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
#endif
		//fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_riemann_analit.dat' using 1:%d w l lw 2 title \"numerical\", 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w l lw 2 title \"exact\"\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);

		fprintf(plot, "\n\n");

	}
	fclose(plot);
	system(name);
}

void gnuplot_analitical_riemann2(int numcells, int* n_r, int* n_u, int* n_p)
{
	FILE *plot;
	char name1[255];
	double dx = LENGTH / double(numcells);
	double time_control[N_smooth];
	double k_step = time_max_array[PROBLEM] / N_smooth;
	printf("step %lf\n", k_step);
	for (int i = 0; i < N_smooth; i++)
		time_control[i] = (i + 1)*k_step;


	for (int i = 0; i < 4; i++) //итерации по характеристикам газа
	{
		if (i == 3) continue;
		/*----------------------------------------*/
		sprintf(name1, "N%03d_P%1d_riemann_analit.plt", numcells, PROBLEM);
		plot = fopen(name1, "w");

		fprintf(plot, "reset\n\
					   clear\n\n\
					   ###################################\n\
					   set term png font \"Times-Roman, 16\"\n\n\
						##################################\n\n");

		for (int k = 0; k < N_smooth; k++)
		{

#ifndef NC
			fprintf(plot, "set xrange[0:1]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
#else
#ifdef NC2
			fprintf(plot, "set xlabel \"(x-x0-Dt)/h\"\n");
			fprintf(plot, "set xrange[-15:15]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
#else
			fprintf(plot, "set xlabel \"x-x0-Dt\"\n");
			fprintf(plot, "set xrange[-0.5:0.5]\n\n"); // забиваем в файл plot все виды функций давления, плотности и скорости в волне разрежения
#endif		
#endif

#ifndef NC
#if PROBLEM!=9
			if (i == 1) fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_SW[i], right_SW[i]);
			else fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_SW[i], right_RW[i]);
#else
			fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_sodd[i], right_sodd[i]);
#endif
#else
			if (PROBLEM == 2) fprintf(plot, "set yrange[%3.2f:%3.2f]\n", left_SW[i], right_SW[i]);
#endif

			if (i == 0) fprintf(plot, "set ylabel \"density\"\n");
			if (i == 1)	fprintf(plot, "set ylabel \"velocity\"\n");
			if (i == 2)	fprintf(plot, "set ylabel \"pressure\"\n");
			if (i == 4) fprintf(plot, "set ylabel \"entropy\"\n");
#ifdef RP
			fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_riemann_analit.dat' using 1:%d w l lw 2 title \"exact\", 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w linespoints pt 4 lt rgb 'blue' title \"num\", 'RP_N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w linespoints pt 7 lt rgb 'green' title \"num_rp\" \n\n",
				numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
#else

#ifndef NC
			if (i != 4)
			{
				fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/N%03d_P%1d_%c_analit_riemann_%6.4lf.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, prop[i], time_control[k]);
				fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_analit_%6.4lf.dat' using 1:%d w l lw 2 title \"exact\", 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%c_%6.4lf.dat' using 1:%d w linespoints pt 7 title \"numerical\"\n\n", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 2, numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, TYPE, time_control[k], i + 2);
			}
			else
			{
				fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/N%03d_P%1d_%c_analit_riemann_%6.4lf.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, prop[i], time_control[k]);
				fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4lf.dat' using 1:%d w linespoints pt 7 title \"numerical\"\n\n", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 2);

				fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/AV_N%03d_P%1d_%c_analit_riemann_%6.4lf.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, prop[i], time_control[k]);
				fprintf(plot, "plot 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4lf.dat' using 1:%d w linespoints pt 7 title \"numerical\"\n\n", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 3);
			}
#else
#ifdef NC2
#ifdef SW_POINTS_PRINT
			if (i == 0) fprintf(plot, "\nset label \"N_{sw} = %d\" at 7, 1.2\n\n", n_r[k]);
			if (i == 1) fprintf(plot, "\nset label \"N_{sw} = %d\" at 7, 0.25\n\n", n_u[k]);
			if (i == 2) fprintf(plot, "\nset label \"N_{sw} = %d\" at 7, 1.3\n\n", n_p[k]);
#endif
#ifdef SW_POINTS_PRINT
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/pr_ND_N%03d_P%1d_%c_analit_riemann_%6.4lf.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, dip[i], time_control[k]);
#else
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/ND_N%03d_P%1d_%c_analit_riemann_%6.4lf.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, dip[i], time_control[k]);
#endif
#else
			fprintf(plot, "set output 'workspace/%03d/%c/P_%1d/NC_N%03d_P%1d_%c_analit_riemann_%6.4lf.png'\n", numcells, prop[i], PROBLEM, numcells, PROBLEM, dip[i], time_control[k]);
#endif
			fprintf(plot, "plot 'workspace/%03d/NC_N%03d_P%1d_SLV%1d_TERM%.0lf_analit_%6.4lf.dat' using 1:%d w l lw 2 title \"exact\", 'workspace/%03d/N%03d_P%1d_SLV%1d_TERM%.0lf_%6.4lf_NC.dat' using 1:%d w linespoints pt 7 title \"numerical\"\n\n", numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 2, numcells, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, time_control[k], i + 2);
#endif


#endif
		}
		//fprintf(plot, "plot 'N%04d_P%1d_SLV%1d_TERM%.0lf_riemann_analit.dat' using 1:%d w l lw 2 title \"numerical\", 'N%04d_P%1d_SLV%1d_TERM%.0lf.dat' using 1:%d w l lw 2 title \"exact\"\n\n", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2, numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM, i + 2);
		fprintf(plot, "\nclear\nreset\n");
		fprintf(plot, "\n\n");

		fclose(plot);
		if (i == 4) system(name1);

	}
}

void output_last_step(int numcells, double dx, double D_analit, double *R, double *U, double *P)
{
	double us, cs, ps, ds, ss, x;
	char FileName[255];
	FILE *fout;
#ifndef RP
	sprintf(FileName, "N%04d_P%1d_SLV%1d_TERM%.0lf.dat", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);
	fout = fopen(FileName, "w");
#else
	sprintf(FileName, "RP_N%04d_P%1d_SLV%1d_TERM%.0lf.dat", numcells, PROBLEM, RUNGE_KUTTA, A_TERM*K_TERM);
	fout = fopen(FileName, "w");
#endif

	for (int i = 0; i < numcells; i++)
	{
		x = i*dx + 0.5*dx;
#ifdef NC
#if (PROBLEM==0)
		x_NC = (i*dx + 0.5*dx) - D_analit*time_max_array[0];  //NC
#elif (PROBLEM==1)
		x_NC = (i*dx + 0.5*dx - 0.1) / (time_max_array[0] * (numb + 1));  //NC
#endif
#endif

		/********************************************************************************************************************************
		| 1 |    2    |   3   |     4    |          5        |         6       |      7      |          8        |      9       |   10  |
		| x | density | speed | pressure | velocity of sound |       entropy   |      term    |	                 |   pg/pg_max  | h_pow |
		********************************************************************************************************************************/
		ds = R[i];
		us = U[i];
		ps = P[i];
		cs = sqrt(GAMMA*P[i] / R[i]);
		ss = log(P[i] / pow(R[i], GAMMA));


#ifndef NC
		fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf \n", x, ds, us, ps, cs, ss);
#else 
		//fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf %30.24lf %30.24lf %30.24lf \n", x_NC, ds, us, ps, cs, es, rp, pg, hh);
		fprintf(fout, "%9.6lf %lf %lf %lf %lf %lf\n", x_layer_NC[i], ds, us, ps, cs, ss);
#endif

	}
}

/**************************************************/