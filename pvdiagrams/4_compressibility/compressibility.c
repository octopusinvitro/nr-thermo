/**********************************************************************************
 * BOYLE: 	Calculates the values of pv/RT and builds the compressibility
			diagram.
 * APPLICATION: Acetone
 **********************************************************************************/	
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../nr/nrutil.c"
#include "../../nr/common.c"

#define R 0.0820336
#define a 15.800111
#define b 0.112363
#define t 10

int main(void) {

	int n, i, j;
	double *v, *p, *pv;
	static double T[t] = {415, 460, 508.1, 530, 548.8, 853.7, 1219.5, 1585.4, 1714.8, 1829.3};

	int numOfLines(char *filename);
	void read1col(int n, double *x, char *filename);

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
	v  = dvector(1,n); 
	p  = dvector(1,n); 
	pv = dvector(1,n); 
	
	/***************************************************
	 * Filling in vectors 
	 ***************************************************/
	read1col(n, v, name);

	for(j=0; j<t; j++) {
		printf(     "\n\nTemperature = %3.3lf K\n", T[j]);
		fprintf(fp, "\n\nTemperature = %3.3lf K\n", T[j]);
		printf(     "\n%10s %10s\n", "p [atm]", "z = pv/RT");
		fprintf(fp, "\n%10s %10s\n", "p [atm]", "z = pv/RT");

		for (i=1; i<=n; i++) {

			p[i]  = (R*T[j] / (v[i]-b)) - (a / pow(v[i],2));
			pv[i] = (p[i] * v[i]) / (R*T[j]);
			printf(     "\n %10.3lf	%10.3lf", p[i], pv[i]);
			fprintf(fp, "\n %10.3lf	%10.3lf", p[i], pv[i]);
		}
	}

	printf(     "\n\n");
	fprintf(fp, "\n\n");

	free_dvector(pv,1,n);
	free_dvector(p,1,n);
	free_dvector(v,1,n);
//	free(T);

	return 0;
}

#undef t
#undef b
#undef a
#undef R
