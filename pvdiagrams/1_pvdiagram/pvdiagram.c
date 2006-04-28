/**********************************************************************************
 * PVDIAGRAM: Program that uses the van der Waals equation to calculate
 *				the curves of the relevant isotherms in the regions that
 *				make sense in a P-V diagra
 * APPLICATION: Plotting p and v in a PV diagram, for the case of acetone.
 **********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../nr/nrutil.c"
#include "../../nr/common.c"

#define R  0.0820336
#define a 15.800111
#define b  0.112363

int main(void) {

	double T;
	double *v, *p;

	int numOfLines(char *filename), N, i;
	void read1col(int n, double *x, char *filename);

	char repeat, name[10];
	FILE *fp;
	if ((fp=fopen("output.txt", "w")) == NULL) {
		printf("\nError opening file\n");
		exit(1);
	}

	/* At an experiment we will have diff v measurements for diff Ts */
	do {
		/***************************************************
		 * Data input
		 ***************************************************/
		printf("\nTemperature (K): ");
		scanf("%lf", &T);
		printf(     "Temperature = %3.3lf K\n", T);
		fprintf(fp, "Temperature = %3.3lf K\n", T);

		printf("\nName of the data file for that T (v in L/mol): ");
		scanf("%s", name);
		N = numOfLines(name);
		printf(     "Number of points = %d\n", N);
		fprintf(fp, "Number of points = %d\n", N);
		getchar();

		/***************************************************
		 * NR vector declaration
		 ***************************************************/
		v = dvector(1,N);
		p = dvector(1,N);

		/***************************************************
		 * Filling in vectors
		 ***************************************************/
		read1col(N, v, name);
		printf(     "\n%15s %15s\n", "v [L/mol]", "p [atm]");
		fprintf(fp, "\n%15s %15s\n", "v [L/mol]", "p [atm]");
		for (i=1; i<=N; i++) {
			p[i] = (R*T / (v[i]-b)) - (a / pow(v[i],2));
			printf(     "\n v[%d] = %3.3lf    p[%d] = %3.3lf", i, v[i], i, p[i]); //%15.6lf
			fprintf(fp, "\n v[%d] = %3.3lf    p[%d] = %3.3lf", i, v[i], i, p[i]);
		}

		// Separator
		printf(     "\n\n");
		fprintf(fp, "\n\n");

		/***************************************************
		* Free memory
		 ***************************************************/
		free_dvector(p,1,N);
		free_dvector(v,1,N);

		printf("Repeat for another temperature? (Y/N): ");
		scanf("%c", &repeat);

	} while(repeat == 'y' || repeat == 'Y');

	getchar();
	return 0;
}

#undef b
#undef a
#undef R
