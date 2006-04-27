/**********************************************************************************
 * SPLINESC: Program that obtains the parameters of a function, using the method
 *           of Cubic splines interpolation.
 *
 * APPLICATION: Resolution of the cross-section according to the Breit-Wigner
 *              formula.
 **********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../nr/nrutil.c"
#include "../../nr/spline.c"
#include "../../nr/splint.c"
#include "../../nr/common.c"

int main(void) {

	int N, n, i, m;
	double *xa, *ya, *x, *y, dx, jmin, jmax, j, yp1, ypn, *y2, gamma;
	
	int numOfLines(char *filename), maximum (double y[], int N);
	void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
	void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
	void read2col(int n, double *x, double *y, char *filename);
	
	char name[10];
	FILE *fp;
	if ((fp=fopen("output.txt", "w")) == NULL) {
		printf("\nError opening file\n");
		exit(1);
	}

	/***************************************************
	 * Data input
	 ***************************************************/
	printf("\nName of data file: ");
	scanf("%s", name);
	N = numOfLines(name);
	printf("N = %d\n", N);
	getchar();

	/***************************************************
	 * NR vector and matrix declaration
	 ***************************************************/
	xa = dvector(1,N);
	ya = dvector(1,N);
	y2 = dvector(1,N);
	
	/***************************************************
	 * Filling in vectors and matrices
	 ***************************************************/
	read2col(N, xa, ya, name);
	printf( "\n\nVectors x and y:\n");
	fprintf(fp, "Vectors x and y:\n");
	for(i=1; i<=N; i++) {
		printf(     "\n E[%d] = %3.3lf    Sig[%d] = %3.3lf", i, xa[i], i, ya[i]);
		fprintf(fp, "\n E[%d] = %3.3lf    Sig[%d] = %3.3lf", i, xa[i], i, ya[i]);
	}

	/***************************************************
	 * Creating the interpolation interval
	 ***************************************************/
	printf("\n\n Lower limit Emin of the interpolation interval: ");
	scanf("%lf", &jmin);	
	printf(" Upper limit Emax of the interpolation interval: ");
	scanf("%lf", &jmax);
	printf(" Step size: ");
	scanf("%lf", &dx);
	
	fprintf(fp, "\n\n Lower limit Emin of the interpolation interval: Emin = %3.3lf", jmin);
	fprintf(fp, "\n Upper limit Emax of the interpolation interval: Emax = %3.3lf", jmax);
	fprintf(fp, "\n Step size: dE = %3.3lf", dx);

	printf(     "\n");
	fprintf(fp, "\n");
	getchar();
		
	n  = (int)( (jmax - jmin) / dx) + 2;	
	x  = dvector(1,n);
	y  = dvector(1,n);
	
	for (i=1, j=jmin; i<=n-1, j<jmax; i++, j=j+dx)
		x[i] = j;
	x[n] = jmax;
	
	printf(     "\nConstructed interval:\n");
	fprintf(fp, "\nConstructed interval:\n");
	for(i=1; i<=n; i++) {
		printf(     "\n e[%d] = %3.3lf", i, x[i]);
		fprintf(fp, "\n e[%d] = %3.3lf", i, x[i]);
	}
	printf(     "\n");
	fprintf(fp, "\n");
	getchar();
	
	/***************************************************
	 * Calculating the first derivative,
	 * evaluating in E1 and E1 and
	 * calling the functions
	 ***************************************************/
	yp1 = -0,00001;
	ypn = -0,5;
	
	spline(xa, ya, N, yp1, ypn, y2);
	for (i=1; i<=n; i++)
		splint(xa, ya, y2, N, x[i], &y[i]);

	printf(     "\nExtrapolated sigma:\n");
	fprintf(fp, "\nExtrapolated sigma:\n");
	for(i=1; i<=n; i++) {
		printf(     "\n sig[%d] = %3.3lf", i, y[i]);
		fprintf(fp, "\n sig[%d] = %3.3lf", i, y[i]);
	}
	printf(     "\n");
	fprintf(fp, "\n");
	getchar();
	
	m      = maximum (y, n);
	gamma  = 2.0*x[m] / sqrt( (y[m]/ya[1]) - 1.0);
	
	printf("\n\tBREIT-WIGNER PARAMETERS\n");
	printf("\n Coefficient: Sigma zero = %3.3lf", y[m]*pow(gamma,2) / 4.0);
	printf("\n Resonance energy: Er    = %3.3lf", x[m]);
	printf("\n Resonande width: Gamma  = %3.3lf", gamma);
	
	fprintf(fp, "\n\n\tBREIT-WIGNER PARAMETERS\n");
	fprintf(fp, "\n Coefficient: Sigma zero = %3.3lf", y[m]*pow(gamma,2) / 4.0);
	fprintf(fp, "\n Resonance energy: Er    = %3.3lf", x[m]);
	fprintf(fp, "\n Resonande width: Gamma  = %3.3lf", gamma);
	
	/***************************************************
	 * Free memory
	 ***************************************************/
	free_dvector(y,  1, n);
	free_dvector(x,  1, n);
	free_dvector(y2, 1, N);
	free_dvector(ya, 1, N);
	free_dvector(xa, 1, N);

	printf("\n\n***************************************************");
	printf("\n\n               END OF CALCULATIONS");
	printf("\n\n***************************************************\n\n");
	return 0;
}
