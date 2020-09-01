#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define Ka 1.5418380 /*wavelength*/
#define rad (M_PI/180)
#define N 4 /*coefficients for analytical approximation to the scattering factors*/
#define M  240
#define Mg 184
#define Al 24
#define Y  32
int main(void) {
	double a = 11.16, b = 11.16, c = 51.62; /*Lattice parameters*/
	double alpha = 90 * rad, beta = 90 * rad, gamma = 120 * rad; /*Lattice parameters*/

																 /*coefficients for analytical approximation to the scattering factors*/
	double
		Mg_a[N] = { 3.49880, 3.83780, 1.32840, 0.849700 },
		Mg_b[N] = { 2.161760, 4.75420, 0.18500, 10.1411 },
		Al_a[N] = { 4.17448, 3.38760, 1.20296, 0.528137 },
		Al_b[N] = { 1.93816, 4.14553, 0.2228753, 8.28524 },
		Zn_a[N] = { 11.9719, 7.38620, 6.46680, 1.39400 },
		Zn_b[N] = { 2.99460, 0.203100, 7.08260, 18.0995 },
		 Y_a[N] = { 17.9268, 9.15310, 1.76795, -33.108 },
		 Y_b[N] = { 1.35417, 11.2145, 2.6599, -0.01319 };

	int i;
	printf("Mg_a[i],""Mg_b[i],""Al_a[i],""Al_b[i],""Zn_a[i],""Zn_b[i],""Y_a[i],""Y_b[i],\n");
	for (i = 0; i < N; i++) {
		printf("%f,""%f,""%f,""%f,""%f,""%f,""%f,""%f,\n"
			, Mg_a[i], Mg_b[i], Al_a[i], Al_b[i], Zn_a[i], Zn_b[i], Y_a[i], Y_b[i]);
	}

	


	int h, k, l, P; /*Miller index*/
	double theta, x, degTheta; /*diffraction angle, args of arcsin*/
	double S, f_Mg, f_Al, f_Zn, f_Y; /*S=sin_theat / lamdba */

									 /*scattering factors*/
	double fs_Mg = 0, f1_Mg = 0.165, f2_Mg = 0.177, C_Mg = 0.858400;
	double fs_Al = 0, f1_Al = 0.204, f2_Al = 0.246, C_Al = 1.115100;
	double fs_Zn = 0, f1_Zn = -0.1612, f2_Zn = 0.678, C_Zn = 1.304100;
	double fs_Y = 0, f1_Y = -0.386, f2_Y = 2.025, C_Y = 1.912130;
	double F_Mg = 0.0, F_Al = 0.0, F_Zn = 0.0, F_Y = 0.0;
	double L, I;


	printf("P,""h,""k," "l,""Rad_theta," "deg_2Theta,""L,""S,"
		"f(s)_Mg," "f(s)_Al," "f(s)_Zn," "f(s)_Y," "f_Mg," "f_Al,"" f_Zn,"" f_Y,"
		" F_Mg," "F_Al," "F_Zn," "F_Y," "Re," "Im," "I,\n");
		


		int e;
		double x1[Mg], y1[Mg], z1[Mg];
	
		int f;
		double x2[Al], y2[Al], z2[Al];

		int g;
		double x3[Y], y3[Y], z3[Y];

	for (h =0; h <= 30; h++) {
		for (k =0; k <= h; k++) {
			for (l = 0; l <= 30; l++) {

				if (h == 0 && k == 00) {
					P = 2;
				}
				else if ((h == k && l == 0) || (k == 0 && l == 0) || (h == 0 && l == 0)) {
					P = 6;
				}
				else if ((h == k && l != 0) || (h == k && l == 0) || (h == 0 && k != l) || (k == 0 && h != l) || (h != k && l == 0) || (h == l && k == 0) || (k == l && h == 0)) {
					P = 12;
				}
				else {
					P = 24;

				}

				x = ((Ka / 2.0) *sqrt((4.0 / 3.0 * ((pow(h, 2.0) + h*k + pow(k, 2.0)) / pow(a, 2.0)) + pow(l, 2.0) / pow(c, 2.0))));

				if (x <= 1 && x >= -1) {

					theta = asin(x);
					degTheta = 2.0 * theta * 180.0 / M_PI;
					L = (1.0 + pow(cos(2.0 * theta), 2.0)) / (2.0 * pow(sin(theta), 2.0) * cos(theta));
					S = pow((sin(theta) / Ka), 2.0);

					/*

					***–³Œø‰»***
					long d;
					for (d = 0; d < N; d++) {
					fs_Mg += ((Mg_a[d]) / exp(Mg_b[d] * S));
					fs_Al += ((Al_a[d]) / exp(Al_b[d] * S));
					fs_Zn += ((Zn_a[d]) / exp(Zn_b[d] * S));
					fs_Y  += (( Y_a[d]) / exp( Y_b[d] * S));
					}

					*/

					fs_Mg = ((Mg_a[0]) / (exp(Mg_b[0] * S))) + ((Mg_a[1]) / (exp(Mg_b[1] * S))) + ((Mg_a[2]) / (exp(Mg_b[2] * S))) + ((Mg_a[3]) / (exp(Mg_b[3] * S))) + C_Mg;
					fs_Al = ((Al_a[0]) / (exp(Al_b[0] * S))) + ((Al_a[1]) / (exp(Al_b[1] * S))) + ((Al_a[2]) / (exp(Al_b[2] * S))) + ((Al_a[3]) / (exp(Al_b[3] * S))) + C_Al;
					fs_Zn = ((Zn_a[0]) / (exp(Zn_b[0] * S))) + ((Zn_a[1]) / (exp(Zn_b[1] * S))) + ((Zn_a[2]) / (exp(Zn_b[2] * S))) + ((Zn_a[3]) / (exp(Zn_b[3] * S))) + C_Zn;
					fs_Y  = (( Y_a[0]) / (exp( Y_b[0] * S))) + (( Y_a[1]) / (exp( Y_b[1] * S))) + (( Y_a[2]) / (exp( Y_b[2] * S))) + (( Y_a[3]) / (exp( Y_b[3] * S))) + C_Y ;

					f_Mg = sqrt(pow((fs_Mg + f1_Mg), 2) + pow(f2_Mg, 2));
					f_Al = sqrt(pow((fs_Al + f1_Al), 2) + pow(f2_Al, 2));
					f_Zn = sqrt(pow((fs_Zn + f1_Zn), 2) + pow(f2_Zn, 2));
					f_Y  = sqrt(pow(( fs_Y +  f1_Y), 2) + pow( f2_Y, 2));

					double Re_Zn = 0.0, Im_Zn = 0.0;
					
			
					
					FILE *datafile1;
					char *filename1;
					filename1 = "AP_Mg.dat"; /* Atomic positions */
				
					double Re_Mg = 0.0, Im_Mg = 0.0;

					datafile1 = fopen(filename1, "r");
						for (e = 0; e < Mg; e++) {
							fscanf(datafile1, "%lf %lf %lf", &x1[e], &y1[e], &z1[e]);

							Re_Mg += cos(2 * M_PI*(h*x1[e] + k*y1[e] + l*z1[e]));
							Im_Mg += sin(2 * M_PI*(h*x1[e] + k*y1[e] + l*z1[e]));
						}

					fclose(datafile1);

					FILE *datafile2;
					char *filename2;
					filename2 = "AP_Al.dat"; /* Atomic positions */

					double Re_Al = 0.0, Im_Al = 0.0;				

					datafile2 = fopen(filename2, "r");
						for (f = 0; f < Al; f++) {
							fscanf(datafile2, "%lf %lf %lf", &x2[f], &y2[f], &z2[f]);
							Re_Al += cos(2 * M_PI*(h*x2[f] + k*y2[f] + l*z2[f]));
							Im_Al += sin(2 * M_PI*(h*x2[f] + k*y2[f] + l*z2[f]));
						}
					fclose(datafile2);


					FILE *datafile3;
					char *filename3;
					filename3 = "AP_Y.dat"; /* Atomic positions */

					double Re_Y = 0.0, Im_Y = 0.0;

					datafile3 = fopen(filename3, "r");
						for (g = 0; g < Y; g++) {
							fscanf(datafile3, "%lf %lf %lf", &x3[g], &y3[g], &z3[g]);
							Re_Y += cos(2 * M_PI*(h*x3[g] + k*y3[g] + l*z3[g]));
							Im_Y += sin(2 * M_PI*(h*x3[g] + k*y3[g] + l*z3[g]));
					
						}
					fclose(datafile3);
					
				

					F_Mg = pow((f_Mg * Re_Mg), 2) + pow((f_Mg * Im_Mg), 2);
					F_Al = pow((f_Al * Re_Al), 2) + pow((f_Al * Im_Al), 2);
					F_Zn = pow((f_Zn * Re_Zn), 2) + pow((f_Zn * Im_Zn), 2);
					F_Y = pow((f_Y  * Re_Y), 2) + pow((f_Y  * Im_Y), 2);

					I = L*(F_Mg + F_Al + F_Y)*S;

					if (I > 1) {
						printf("%d,%d,%d,%d,", P, h, k, l);
						printf("%7.5f, %7.5f, %f, %f, ", theta, degTheta, L, S);
						printf("%7.5f, %7.5f, %7.5f, %7.5f,", fs_Mg, fs_Al, fs_Zn, fs_Y);
						printf("%7.5f, %7.5f, %7.5f, %7.5f,", f_Mg, f_Al, f_Zn, f_Y);
						printf("%7.5f," "%7.5f," "%7.5f," "%7.5f,", F_Mg, F_Al, F_Zn, F_Y);
						printf("%f,""%f,""%f \n", Re_Mg, Im_Mg, I);
						
					
					}
	

				}


			}

		}
	}
						printf("\n");


	for (e = 0; e < 184; e++) {
		printf("Mg[%d] %f %f %f\n", e, x1[e], y1[e], z1[e]);
	}
		

	for (f = 0; f<24; f++)
		printf("Al[%d] %f %f %f\n",f, x2[f], y2[f], z2[f]);


	for (g = 0; g<32; g++)
		printf("Y[%d] %f %f %f\n",g, x3[g], y3[g], z3[g]);
	
	return 0;
}
