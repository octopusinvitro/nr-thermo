/**********************************************************************************
 * GIBBS POTENTIAL:	Calculates the excess Gibbs potential from vapor-liquid
					equilibrium data.
 * APPLICATION:	Study of the system 1,3-dioxolane + hexane at 308.15 K.
 **********************************************************************************/	
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../nr/nrutil.c"
#include "../../nr/common.c"

#define R 8.31451
#define T 308.15

int main(void) {

	int n, i;
	double *x, *y, *p;
	double  v1o = 71, 
			v2o = 132, 
			p1o = 21.457, 
			p2o = 30.636, 
			B11 = -1439, 
			B22 = -1726, 
			B12 = -1396;
	double *m1E, *m2E, *gE, d12 = 2*B12 - B11 - B22;
	
	char name[10];
	FILE *fp;
	if ((fp=fopen("output.txt", "w")) == NULL) {
		printf("\nError opening file\n");
		exit(1);
	}

	/***************************************************
	 * Data input
	 ***************************************************/
	printf("\nName of the data file (v in L/mol): ");
	scanf("%s", name);
	n = numOfLines(name);
	printf(		"Number of points = %d", n);
	fprintf(fp, "Number of points = %d", n);

	/***************************************************
	 * NR vector declaration
	 ***************************************************/
	x   = dvector(1,n); 
	y   = dvector(1,n); 
	p   = dvector(1,n); 
	m1E = dvector(1,n);
	m2E = dvector(1,n);
	gE  = dvector(1,n); 
	
	/***************************************************
	 * Filling in vectors 
	 ***************************************************/
	read3col(n, x, y, p, name);

	printf(     "\n\nExperimental values: molar fractions (V and L phases), vapor pressure:\n");
	fprintf(fp, "\n\nExperimental values: molar fractions (V and L phases), vapor pressure:\n");
	printf(     "\n%10s %10s %10s\n", "x1", "y1", "P [KPa]");	
	fprintf(fp, "\n%10s %10s %10s\n", "x1", "y1", "P [KPa]");	

	for(i=1; i<=n; i++) {
		printf(     "\n%10.3lf %10.3lf %10.3lf", x[i], y[i], p[i]);
		fprintf(fp, "\n%10.3lf %10.3lf %10.3lf", x[i], y[i], p[i]);
	}
		
	for (i=1; i<=n; i++) {

		if (x[i] == 0)
			m1E[i] = 0;
		else
			m1E[i] = R*T*log((y[i]*p[i]) / (x[i]*p1o))
					 + ((p[i] - p1o)*(B11 - v1o)
					 + d12*pow(y[i], 2)*p[i])*0.001;

		if (x[i] == 1)
			m2E[i] = 0;
		else
			m2E[i] = R*T*log(((1 - y[i])*p[i]) / ((1 - x[i])*p2o))
					 + ((p[i] - p2o)*(B22 - v2o)
					 + d12*pow(1 - y[i], 2)*p[i])*0.001;

		gE[i]  = x[i]*m1E[i] + (1 - x[i])*m2E[i];
	}
	
	printf(     "\n\nChemical and excess Gibbs potentials:\n");
	fprintf(fp, "\n\nChemical and excess Gibbs potentials:\n");
	printf(     "\n%10s %10s %10s %10s\n", "x1", "m1E [J/mol]", "m2E [J/mol]", "gE [J/mol]");	
	fprintf(fp, "\n%10s %10s %10s %10s\n", "x1", "m1E [J/mol]", "m2E [J/mol]", "gE [J/mol]");

	for(i=1; i<=n; i++) {
		printf(     "\n %10.3lf %10.3lf %10.3lf %10.3lf", x[i], m1E[i], m2E[i], gE[i]);
		fprintf(fp, "\n %10.3lf %10.3lf %10.3lf %10.3lf", x[i], m1E[i], m2E[i], gE[i]);
	}
	
	printf("\n\n");
	free_dvector(gE,1,n);
	free_dvector(m2E,1,n);
	free_dvector(m1E,1,n);
	free_dvector(p,1,n);
	free_dvector(y,1,n);
	free_dvector(x,1,n);
	return 0;
}

#undef T
#undef R
