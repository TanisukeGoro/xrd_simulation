#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define Ka 1.541838                                         /*wavelength*/
#define rad (M_PI/180)
#define N 10


int main(int argc, char *argv[]) {

	double z[N], x[N], y[N], a = 0.0, c = 0.0;
	double xx = 0.0, xy = 0.0, yy = 0.0, xz = 0.0, yz = 0.0;
	double theta[N]; //Diffraction angle
	int h[N], k[N], l[N];
	int i, m, b, q;
	double resum=0.0, stdev=0.0; //residual sum of squares, Saple standard deviation
	double Ga=0.0, Gc=0.0; //Matrix G^-1 _11, _22
	double Sighat=0.0; //sigma hat
	double t=2.306; //Student's t-distribution_95%ture
	double Cla=0.0, Clc=0.0; //Confidence interval
	
	printf("[Wave length]\n %9.7f\n", Ka);
	
	FILE *datafile;
	char *filename;
	filename = "Mg75Al10Y15.dat"; // hkl, diffraction angle

	datafile = fopen(filename, "r");

	if ((datafile = fopen(filename, "r")) == NULL) {
		printf("Can not open file - %s\n", filename);
		return (-1);
	}
	else {

		for (m = 0; m < N; m++) {

			fscanf(datafile, "%d %d %d %lf", &h[m], &k[m], &l[m], &theta[m]);

		}

		printf("file scan complete!!\n");
		printf("[%s]\n h, l, l, 2theta\n", filename);

	}
	fclose(datafile);

	for (q = 0; q < N; q++) {
		printf("%d, %d, %d, %7.5f\n", h[q], k[q], l[q], theta[q]);

	}



	for (b = 0; b < N; b++) {

		x[b] = (pow(Ka, 2) / 3.0)*(pow(h[b], 2) + h[b] * k[b] + pow(k[b], 2));
		y[b] = (pow(Ka, 2) / 4.0)*pow(l[b], 2);
		z[b] = pow((sin(rad*theta[b] / 2.0)), 2);

	}
	printf("x[],y[],z[]\n");
	for (b = 0; b < N; b++) {
		printf("%f, %f, %f\n", x[b], y[b], z[b]);
	}


		for (i = 0; i < N; i++) {
			xz = xz + x[i] * z[i];
			yz = yz + y[i] * z[i];
			xx = xx + x[i] * x[i];
			yy = yy + y[i] * y[i];
			xy = xy + x[i] * y[i];

		}

	//printf("xz, yz, xx, yy, xy\n");  ****debug****
	//printf("%7.5f, %7.5f, %7.5f, %7.5f, %7.5f\n", xz, yz, xx, yy, xy); ****debug****

	a = sqrt((xx*yy - pow(xy, 2.0)) / (yy*xz - xy*yz));
	c = sqrt((xx*yy - pow(xy, 2.0)) / (xx*yz - xy*xz));
	
		for(i=0; i<N; i++){
	resum=resum + powl((z[i] - (a * x[i] + c * y[i] )), 2 );
		}

	Sighat= sqrt(resum/(N-2.0));
	Ga=sqrt(yy/(xx*yy*pow(xy,2.0)));
	Gc=sqrt(xx/(xx*yy*pow(xy,2.0)));
	Cla=t*Ga*Sighat/sqrt(N);	
	Clc=t*Gc*Sighat/sqrt(N);
	

	printf("[Lattic Parameters]\n %f, %f\n", a, c);
	printf("[Confidence interval]\n %f, %f\n", Cla, Clc);

	return 0;

}