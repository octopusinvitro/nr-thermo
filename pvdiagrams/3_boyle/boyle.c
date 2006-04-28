/**********************************************************************************
 * BOYLE: 	Calculates the curve that connects the minima of horizontal slope 
			of the pv-p curve.
 * APPLICATION: Acetone
 **********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../nr/nrutil.h"
#include "../../nr/nrutil.c"
#include "../../nr/zrhqr.c"
#include "../../nr/balanc.c"
#include "../../nr/hqr.c"
#include "../../nr/common.c"

#define R 0.0820336
#define a 15.800111
#define b 0.112363
#define t 10
#define M 2	 
#define N (M+1) 

int main(void) {

	int n, i, j, l;
	double *v, *p, *pv, *pvmin, *pmin, lower, c[N], *rti, *rtr, *rt, *pboyle;
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
	v      = dvector(1,n);
	p      = dvector(1,n);
	pv     = dvector(1,n);
	pvmin  = dvector(1,t);
	pmin   = dvector(1,t);
	rti    = dvector(1,M);
	rtr    = dvector(1,M);
	rt     = dvector(1,2*n);
	pboyle = dvector(1,n);
	
	/***************************************************
	 * Filling in vectors 
	 ***************************************************/
	read1col(n, v, name);

	for(j=0; j<t; j++) {

		l     = 0;
		lower = 1000000000;
	
		printf(     "\n\nTemperature = %3.3lf K\n", T[j]);
		fprintf(fp, "\n\nTemperature = %3.3lf K\n", T[j]);
		printf(     "\n%15s %15s %15s\n", "v [L/mol]", "p [atm]", "pv [atmL/mol]");
		fprintf(fp, "\n%15s %15s %15s\n", "v [L/mol]", "p [atm]", "pv [atmL/mol]");

		for (i=1; i<=n; i++) {

			p[i]  = (R*T[j] / (v[i]-b)) - (a / pow(v[i],2));
			pv[i] = p[i] * v[i];
			printf(     "\n %10.3lf	%10.3lf	%10.3lf", v[i], p[i], pv[i]);
			fprintf(fp, "\n %10.3lf	%10.3lf	%10.3lf", v[i], p[i], pv[i]);

			c[0] = 2*a*b*p[i];
			c[1] = -a;
			c[2] = b;
			zrhqr(c, M, rtr, rti);

			if ((rti[1] == 0.0) && 
				(rti[2] == 0.0)) {
				if ((rtr[1] != 0.0) && 
					(rtr[2] != 0.0)) {
					l++;
					rt[l]     = rtr[1];
					rt[l+n]   = rtr[2];
					pboyle[l] = p[i];
				}
			}

			if(pv[i] < lower) {
				lower = pv[i];
				pmin[j+1] = p[i];
			}
		}

		pvmin[j+1] = lower;
				
		printf(     "\n\n Roots of the Boyle curve:\n");
		fprintf(fp, "\n\n Roots of the Boyle curve:\n");
		printf(     "\n%15s %15s %15s\n", "p [atm]", "pv1 [atm/mol]", "pv2 [atmL/mol]");
		fprintf(fp, "\n%15s %15s %15s\n", "p [atm]", "pv1 [atm/mol]", "pv2 [atmL/mol]");

		for (i=1; i<=l; i++) {
			printf(     "\n %10.3lf	%10.3lf	%10.3lf", pboyle[i], rt[i], rt[i+n]);
			fprintf(fp, "\n %10.3lf	%10.3lf	%10.3lf", pboyle[i], rt[i], rt[i+n]);
		}

		// Separator
		printf(     "\n\n------------------------------------");
		fprintf(fp, "\n\n------------------------------------");
		getchar();
	}
		
	printf(     "\nMinima of the pvp-p curves\n");
	fprintf(fp, "\nMinima of the pvp-p curves\n");
	printf(     "\n%15s %10s\n", "pmin [atm]", "pvmin [atmL/mol]");
	fprintf(fp, "\n%15s %10s\n", "pmin [atm]", "pvmin [atmL/mol]");

	for (i=1; i<=t; i++) {	
		printf(     "\n %10.3lf	%10.3lf", pmin[i], pvmin[i]);
		fprintf(fp, "\n %10.3lf	%10.3lf", pmin[i], pvmin[i]);
	}

	lower = 0.0;	/*	p = 0, (to avoid declaring new variable) */
	c[0] = 2*a*b*lower;
	c[1] = -a;
	c[2] = b;
	zrhqr(c, M, rtr, rti);

	for (i=1; i<=M; i++) {
		if ((rti[i] == 0.0) && 
			(rtr[i] != 0.0)) {	/* Just in case */
			printf(     "\n\n Calculated Boyle temperature            = %3.3lf", rtr[i] / R);
			fprintf(fp, "\n\n Calculated Boyle temperature            = %3.3lf", rtr[i] / R);
		}
	}

	printf(     "\n Boyle's limit temperature (a/bR)        = %3.3lf", a / (b*R));
	fprintf(fp, "\n Boyle's limit temperature (a/bR)        = %3.3lf", a / (b*R));
	printf(     "\n Boyle's temperature from Tc (27Tc/8) up = %3.3lf", (27.0*508.1) / 8.0);
	fprintf(fp, "\n Boyle's temperature from Tc (27Tc/8) up = %3.3lf", (27.0*508.1) / 8.0);
	printf(     "\n\n");
	fprintf(fp, "\n\n");

	free_dvector(pboyle,1,n);	
	free_dvector(rt,1,2*n);	
	free_dvector(rtr,1,M);
	free_dvector(rti,1,M); 
	free_dvector(pmin,1,t);		
	free_dvector(pvmin,1,t);	
	free_dvector(pv,1,n);
	free_dvector(p,1,n);
	free_dvector(v,1,n);
//	free(T);
//	free(c);

	return 0;
}

#undef N
#undef M
#undef t
#undef b
#undef a
#undef R
