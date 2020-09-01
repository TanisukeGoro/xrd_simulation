

int input_isbn(void);

int main() 
{
	system("XRDTi.txt");
	Sleep(1000);
	input_isbn();
	getch();
	printf("終わりだよ！");
	return 0;
}

int input_isbn(void)
{
	char isbn[64];

	printf("ISBNを入力してください。　\n");
		scanf("%s", isbn); //文字列
	printf("入力された番号は、%sです。\n", isbn);
	return 0;

}
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include<conio.h>
#include <windows.h>

#define Ka 1.540562 /*wavelength*/
#define rad (M_PI/180)
#define N 4 /*coefficients for analytical approximation to the scattering factors*/
#define M  240
#define Mg 184
#define Al 24
#define Y  32
#define L12 4
//Declaration of Atomic Species
char *Mg_ATOMIC_SPECIES = "Mg";//Mg or Mg+2
char *Al_ATOMIC_SPECIES = "Al+3";	 //Al or Al+3
char *Zn_ATOMIC_SPECIES = "Zn+2";  //Zn or Zn+2
char *Y_ATOMIC_SPECIES = "Y+3";  //Y  or Y+3

								 //Atomic Species Data
double
a_1[N] = { 5.42040, 2.17350, 1.22690, 2.370730 },//Mg
b_1[N] = { 2.82750, 79.2611, 0.380800, 7.19370 },//Mg
c_1 = 0.858400,						   //Mg
A_1[N] = { 3.49880, 3.83780, 1.32840, 0.849700 },//Mg+2
B_1[N] = { 2.161760, 4.75420, 0.18500, 10.1411 },//Mg+2
C_1 = 0.485300,								   //Mg+2

a_2[N] = { 6.42020, 1.90020, 1.59360, 1.96460 }, //Al
b_2[N] = { 3.03870, 0.742600, 31.5472, 85.0886 },//Al
c_2 = 1.115100,									 //Al
A_2[N] = { 4.17448, 3.38760, 1.20296, 0.528137 },//Al+3
B_2[N] = { 1.93816, 4.14553, 0.228753, 8.28524 },//Al+3		
C_2 = 0.706786,									 //Al+3

a_3[N] = { 14.0743, 7.03180, 5.16520, 2.41000 },//Zn
b_3[N] = { 3.26550,0.233300, 10.3163, 58.7097 },//Zn
c_3 = 1.30410,									//Zn
A_3[N] = { 11.9719, 7.38620, 6.46680, 1.39400 },//Zn+2
B_3[N] = { 2.99460, 0.20310, 7.08260, 18.0995 },//Zn+2
C_3 = 0.780700,									//Zn+2

a_4[N] = { 17.7760, 10.2946, 5.72629, 3.26588 },//Y
b_4[N] = { 1.40290, 12.8006, 0.125599, 104.354 },//Y
c_4 = 1.912130,									//Y
A_4[N] = { 17.9268, 9.15310, 1.76795, -33.108 },//Y+3
B_4[N] = { 1.35417, 11.2145, 22.6599, -0.01319 },//Y+3
C_4 = 40.2602;									//Y+3

int main(int argc, char *argv[]) {
	
	
	
	double a = 11.262, b = 11.262, c = 52.544; /*Lattice parameters*/
	double alpha = 90 * rad, beta = 90 * rad, gamma = 120 * rad; /*Lattice parameters*/

																 //ATOMIC_SPECIES
	char *st_Mg = "Mg";
	char *st_Al = "Al";
	char *st_Zn = "Zn";
	char *st_Y = "Y";
	int k_1, k_2, k_3, k_4;
	k_1 = strcmp(Mg_ATOMIC_SPECIES, st_Mg);
	k_2 = strcmp(Al_ATOMIC_SPECIES, st_Al);
	k_3 = strcmp(Zn_ATOMIC_SPECIES, st_Zn);
	k_4 = strcmp(Y_ATOMIC_SPECIES, st_Y);
	double Mg_a[N], Mg_b[N], C_Mg;
	double Al_a[N], Al_b[N], C_Al;
	double Zn_a[N], Zn_b[N], C_Zn;
	double Y_a[N], Y_b[N], C_Y;


	if (k_1 == 0) {
		memcpy(Mg_a, a_1, sizeof(a_1));
		memcpy(Mg_b, b_1, sizeof(b_1));
		C_Mg = c_1;
	}
	else {
		memcpy(Mg_a, A_1, sizeof(A_1));
		memcpy(Mg_b, B_1, sizeof(B_1));
		C_Mg = C_1;
	}
	if (k_2 == 0) {
		memcpy(Al_a, a_2, sizeof(a_2));
		memcpy(Al_b, b_2, sizeof(b_2));
		C_Al = c_2;
	}
	else
	{
		memcpy(Al_a, A_2, sizeof(A_2));
		memcpy(Al_b, B_2, sizeof(B_2));
		C_Al = C_2;
	}

	if (k_3 == 0) {
		memcpy(Zn_a, a_3, sizeof(a_3));
		memcpy(Zn_b, b_3, sizeof(b_3));
		C_Zn = c_3;
	}
	else
	{
		memcpy(Zn_a, A_3, sizeof(A_3));
		memcpy(Zn_b, B_3, sizeof(B_3));
		C_Zn = C_3;
	}

	if (k_4 == 0) {
		memcpy(Y_a, a_4, sizeof(a_4));
		memcpy(Y_b, b_4, sizeof(b_4));
		C_Y = c_4;
	}
	else
	{
		memcpy(Y_a, A_4, sizeof(A_4));
		memcpy(Y_b, B_4, sizeof(B_4));
		C_Y = C_4;
	}

	printf("USED_ATOMIC_SPECIES\n%s,,,%s,,,%s,,,%s\n", Mg_ATOMIC_SPECIES, Al_ATOMIC_SPECIES, Zn_ATOMIC_SPECIES, Y_ATOMIC_SPECIES);
	int i, j = 1;

	for (i = 0; i < N; i++) {

		printf("Mg_a[%d]=,%f,Mg_b[%d]=,%f,Al_a[%d]=,%f,Al_b[%d]=,%f,Zn_a[%d]=,%f,Zn_b[%d]=,%f,Y_a[%d]=,%f,Y_b[%d]=,%f,\n"
			, j, Mg_a[i], j, Mg_b[i], j, Al_a[i], j, Al_b[i], j, Zn_a[i], j, Zn_b[i], j, Y_a[i], j, Y_b[i]);
		++j;
	}
	printf("C_Mg=,%f,,,C_Al=,%f,,,C_Zn=,%f,,,C_Y=,%f\n", C_Mg, C_Al, C_Zn, C_Y);

	int h, k, l, P; /*Miller index*/
	double theta, d, x, degTheta, q; /*diffraction angle, args of arcsin*/
	double S, f_Mg, f_Al, f_Zn, f_Y; /*S=sin_theat / lamdba */

									 /*scattering factors*/
	double fs_Mg = 0, f1_Mg = 0.165, f2_Mg = 0.177;
	double fs_Al = 0, f1_Al = 0.204, f2_Al = 0.246;
	double fs_Zn = 0, f1_Zn = -1.612, f2_Zn = 0.678;
	double fs_Y = 0, f1_Y = -0.386, f2_Y = 2.025;
	double F_Mg = 0.0, F_Al = 0.0, F_Zn = 0.0, F_Y = 0.0;
	double F = 0.0;
	double L, I;
	double occ_Mg = 0.0, occ_Al = 0.0, occ_Y = 0.0;
	/***************************************************************************************/

	FILE *datafile5;
	char *filename5;
	filename5 = "occ_Al_Y.dat"; /* Atomic positions */
	datafile5 = fopen(filename5, "r");
	fscanf(datafile5, "%lf %lf %lf", &occ_Mg, &occ_Al, &occ_Y);

	fclose(datafile5);

	/****************************************************************************************/
	printf("occ_Mg, occ_Al, occ_Y\n");

	printf("%f, %f, %f\n", occ_Mg, occ_Al, occ_Y);

	printf("P,""h,""k," "l,""D-valu,""Rad_theta," "L,""S,"
		"f(s)_Mg," "f(s)_Al," "f(s)_Y," "f_Mg," "f_Al,"" f_Y,"
		" F," "Re_Mg," "Im_Mg," "deg_2Theta,""I,\n");



	int e;
	double x1[Mg], y1[Mg], z1[Mg];

	int f;
	double x2[Al], y2[Al], z2[Al];

	int g;
	double x3[Y], y3[Y], z3[Y];

	int p;
	double x4[L12], y4[L12], z4[L12];

	for (h = -60; h <= 60; h++) {
		for (k = -60; k <= 60; k++) {
			for (l = -60; l <= 60; l++) {

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

				P = 1;//no multi

				x = ((4.0 / 3.0) * ((h*h + h*k + k*k) / (a*a)) + (l*l) / (c*c));

				d = sqrt(1.0 / x);// d-value

				q = (Ka / 2.0) *sqrt(x);

				S = 0.0;
				L = 0.0;

				if (q <= 1 && q >= -1) {

					theta = asin(q);
					degTheta = 2.0 * theta * 180.0 / M_PI;
					//L = (1.0 + pow(cos(2.0 * theta), 2.0)) / (2.0 * pow(sin(theta), 2.0) * cos(theta));
					L = (1.0 + pow(cos(2.0 * theta), 2.0)) / (sin(2 * theta)*sin(theta));
					S = pow((sin(theta) / Ka), 2.0);


					fs_Mg = ((Mg_a[0] / exp(Mg_b[0] * S)) + (Mg_a[1] / exp(Mg_b[1] * S)) + (Mg_a[2] / exp(Mg_b[2] * S)) + (Mg_a[3] / exp(Mg_b[3] * S))) + C_Mg;
					fs_Al = ((Al_a[0] / exp(Al_b[0] * S)) + (Al_a[1] / exp(Al_b[1] * S)) + (Al_a[2] / exp(Al_b[2] * S)) + (Al_a[3] / exp(Al_b[3] * S))) + C_Al;
					fs_Y = ((Y_a[0] / exp(Y_b[0] * S)) + (Y_a[1] / exp(Y_b[1] * S)) + (Y_a[2] / exp(Y_b[2] * S)) + (Y_a[3] / exp(Y_b[3] * S))) + C_Y;

					f_Mg = sqrt(pow((fs_Mg + f1_Mg), 2) + pow(f2_Mg, 2));
					f_Al = sqrt(pow((fs_Al + f1_Al), 2) + pow(f2_Al, 2));
					f_Y = sqrt(pow((fs_Y + f1_Y), 2) + pow(f2_Y, 2));

					/*

					***無効化***
					long d;
					for (d = 0; d < N; d++) {
					fs_Mg += ((Mg_a[d]) / exp(Mg_b[d] * S));
					fs_Al += ((Al_a[d]) / exp(Al_b[d] * S));
					fs_Zn += ((Zn_a[d]) / exp(Zn_b[d] * S));
					fs_Y  += (( Y_a[d]) / exp( Y_b[d] * S));
					}

					*/

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


					FILE *datafile4;
					char *filename4;
					filename4 = "AP_L12.dat"; /* Atomic positions */

					double Re_L12 = 0.0, Im_L12 = 0.0;

					datafile4 = fopen(filename4, "r");
					for (p = 0; p < L12; p++) {
						fscanf(datafile4, "%lf %lf %lf", &x4[p], &y4[p], &z4[p]);
						Re_L12 += cos(2 * M_PI*(h*x4[p] + k*y4[p] + l*z4[p]));
						Im_L12 += sin(2 * M_PI*(h*x4[p] + k*y4[p] + l*z4[p]));

					}
					fclose(datafile4);

					F = pow((f_Mg * Re_Mg + f_Al*Re_Al + f_Y*Re_Y + (f_Mg*occ_Mg + f_Al*occ_Al + f_Y*occ_Y)*Re_L12), 2.0) + pow((f_Mg * Im_Mg + f_Al*Im_Al + f_Y*Im_Y + ((f_Mg*occ_Mg + f_Al*occ_Al + f_Y*occ_Y)*Im_L12)), 2.0);
					//F_Al = powl((f_Al * Re_Al), 2.0) + powl((f_Al * Im_Al), 2.0);
					//F_Zn = powl((f_Zn * Re_Zn), 2.0) + powl((f_Zn * Im_Zn), 2.0);
					//F_Y = powl((f_Y  * Re_Y), 2.0) + powl((f_Y  * Im_Y), 2.0);

					I = L*F;

					if (I > 1) {




						printf("%d,%d,%d,%d,", P, h, k, l);
						printf("%f,%7.5f, %f, %f, ", d, theta, L, S);
						printf("%7.5f, %7.5f, %7.5f,", fs_Mg, fs_Al, fs_Y);
						printf("%7.5f, %7.5f, %7.5f,", f_Mg, f_Al, f_Y);
						printf("%7.5f,", F);
						printf("%f,""%f,""%7.5f,""%f \n", Re_Mg, Im_Mg, degTheta, I);


					}


				}


			}

		}
	}
	printf("\n");

	j = 1;
	for (e = 0; e < 184; e++) {
		printf("Mg[%d] %f %f %f\n", j, x1[e], y1[e], z1[e]);
		j++;
	}

	j = 1;
	for (f = 0; f < 24; f++) {
		printf("Al[%d] %f %f %f\n", j, x2[f], y2[f], z2[f]);
		j++;
	}
	j = 1;
	for (g = 0; g < 32; g++) {
		printf("Y[%d] %f %f %f\n", j, x3[g], y3[g], z3[g]);
		j++;
	}
	j = 1;
	for (g = 0; g < 4; g++) {
		printf("L12[%d] %f %f %f\n", j, x4[g], y4[g], z4[g]);
		j++;
	}


	return 0;
}
