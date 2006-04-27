/**********************************************************************************
 * BINODAL: Calculates the total area contained by the isotherm and p_0 = cte
 * APPLICATION: Acetone
 **********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../../nr/nrutil.h"
#include "../../nr/nrutil.c"
#include "../../nr/zbrac.c"
#include "../../nr/zbrent.c"
#include "../../nr/zrhqr.c"
#include "../../nr/balanc.c"
#include "../../nr/hqr.c"
#include "../../nr/common.c"

#define R  0.0820336
#define a 15.800111
#define b  0.112363
#define M 3	 
#define N (M+1)		
#include "func.c"

int main(void) {

	int n, i, j, k;
	double *v, *p;
	double T, p1, p2, p3, po, coeff[N], *rti, *rtr, aux, sum;

	int numOfLines(char *filename);
	void read1col(int n, double *x, char *filename);

	char repeat, name[10];
	FILE *fp;
	if ((fp=fopen("output.txt", "w")) == NULL) {
		printf("\nError opening file\n");
		exit(1);
	}
		
	do {
		/***************************************************
		 * Data input
		 ***************************************************/
		printf("\nTemperature (K): ");
		scanf("%lf", &T);
		printf(	    "Temperature = %3.3lf K\n", T);
		fprintf(fp, "Temperature = %3.3lf K\n", T);

		printf("\nName of the data file for that T (v in L/mol): ");
		scanf("%s", name);
		n = numOfLines(name);
		printf(	    "Number of points = %d\n", n);
		fprintf(fp, "Number of points = %d\n", n);
		getchar();
		
		/***************************************************
		 * NR vector declaration
		 ***************************************************/
		v	= dvector(1,n);
		p	= dvector(1,n);
		rti = dvector(1,M);
		rtr = dvector(1,M);
		
		/***************************************************
		 * Filling in vectors 
		 ***************************************************/
		read1col(n, v, name);
		printf(     "\n%10s %10s\n", "v [L/mol]", "p [atm]");	
		fprintf(fp, "\n%10s %10s\n", "v [L/mol]", "p [atm]");	
		for (i=1; i<=n; i++) {
			p[i] = (R*T / (v[i]-b)) - (a / pow(v[i],2));
			printf(     "\n %10.3lf	%10.3lf", v[i], p[i]);
			fprintf(fp, "\n %10.3lf	%10.3lf", v[i], p[i]);
		}

		// Separator
		printf(     "\n");
		fprintf(fp, "\n");
	 
		// For acetone:
		if (T >= 508.1) {
			printf(     "\n Ideal region, above critical p -> incompressible gas");
			fprintf(fp, "\n Ideal region, above critical p -> incompressible gas");
			printf(     "\n Critical point: (pc, Tc) = (46.35 atm, 508.1 K)");
			fprintf(fp, "\n Critical point: (pc, Tc) = (46.35 atm, 508.1 K)");
		  
			coeff[0] = -a*b;
			coeff[1] = a;
			coeff[2] = -(46.35*b + R*508.1);
			coeff[3] = 46.35;			
			zrhqr(coeff, M, rtr, rti);
		
			/* Although we should have just 1 root, we lack that precision at the critical point. */
			printf(     "\n Roots of the curve at Tc = 508.1 K:\n");
			fprintf(fp, "\n Roots of the curve at Tc = 508.1 K:\n");
			printf(     "\n%10s %10s\n", "v Real", "v Imaginary");
			fprintf(fp, "\n%10s %10s\n", "v Real", "v Imaginary");
			for (i=1; i<=M; i++) {
				if(rti[i] == 0.0) {
					printf(     "%10.3f %10.3f\n", rtr[i], rti[i]);
					fprintf(fp, "%10.3f %10.3f\n", rtr[i], rti[i]);
				}
			}
		}
		
		/**************************************************************** 
		 * We calculate the roots of dp/dv = 0 (third grade polynomial)
		 ****************************************************************/
		else {
			coeff[0] = -2*a*b*b;
			coeff[1] = 4*a*b;
			coeff[2] = -2*a;
			coeff[3] = R*T;		 
			zrhqr(coeff, M, rtr, rti);

			printf(     "\n%10s %10s\n", "v Real", "v Imaginary");
			fprintf(fp, "\n%10s %10s\n", "v Real", "v Imaginary");
			for (j=1; j<=M; j++) {
				printf(     "%10.3f  %10.3f\n", rtr[j], rti[j]);
				fprintf(fp, "%10.3f  %10.3f\n", rtr[j], rti[j]);
			}
		
			/* Now we calculate p. Sort the 3 roots from lower to higher */
			for (k=1; k<M; k++) {
				for (j=k+1; j<=M; j++) {
					if(rtr[k] > rtr[j]) {
						aux    = rtr[i];
						rtr[k] = rtr[j];
						rtr[k] = aux;
					}
				} 
			}
		  
			p1 = ((R*T) / (rtr[1] - b)) - (a / pow(rtr[1], 2));
			p2 = ((R*T) / (rtr[2] - b)) - (a / pow(rtr[2], 2));
			p3 = ((R*T) / (rtr[3] - b)) - (a / pow(rtr[3], 2));
			printf(     "\n p1 = %lf, p2 = %lf, p3 = %lf\n", p1, p2, p3);
			fprintf(fp, "\n p1 = %lf, p2 = %lf, p3 = %lf\n", p1, p2, p3);
			
			/* Take the two p's evaluated in the smallest v's */
			po = zbrent(func, p1, p2, 0.01, T);
			printf(     "\n p0 = %lf\n", po);
			fprintf(fp, "\n p0 = %lf\n", po);
			
			/* With p0 we calculate the 2 points of the binodal for this T */
			coeff[0] = -a*b;
			coeff[1] = a;
			coeff[2] = -(po*b + R*T);
			coeff[3] = po;			
			zrhqr(coeff, M, rtr, rti);

			printf(     "\n Roots of the curve at T = %lf:\n", T);
			fprintf(fp, "\n Roots of the curve at T = %lf:\n", T);
			printf(     "\n%10s %5s\n", "v Real", "v Imaginary");
			fprintf(fp, "\n%10s %5s\n", "v Real", "v Imaginary");
			for (i=1; i<=M; i++) {
				 if (rti[i] == 0.0) {
					printf(     "%5d %3.3f %3.3f\n", i, rtr[i], rti[i]);
					fprintf(fp, "%5d %3.3f %3.3f\n", i, rtr[i], rti[i]);
				}
			}

			/* We take only the real ones */
			if ((rti[1] == 0.0) && 
				(rti[2] == 0.0) && 
				(rti[3] == 0.0)) {
				for (k=1; k<M; k++) {
					for (j=k+1; j<=M; j++) {
						if(rtr[k] > rtr[j]) {
							 aux    = rtr[i];
							 rtr[k] = rtr[j];
							 rtr[k] = aux;
						}
					}
				}
				
				printf(     "\n\nPoints where the binodal curve touches p0:\n");
				fprintf(fp, "\n\nPoints where the binodal curve touches p0:\n");
				printf(     "\np0 = %lf atm\t v1 = %lf L/mol\t v2 = %lf L/mol", po, rtr[1], rtr[3]);
				fprintf(fp, "\np0 = %lf atm\t v1 = %lf L/mol\t v2 = %lf L/mol", po, rtr[1], rtr[3]);
				printf(     "\n\n");
				fprintf(fp, "\n\n");
			}
		
			else {
				printf(     "\nThe roots are complex.\n");
				fprintf(fp, "\nThe roots are complex.\n");
			}
		}

		// Separator
		printf(     "\n------------------------------------\n\n");
		fprintf(fp, "\n------------------------------------\n\n");

	
		/***************************************************
	 	* Free memory
		 ***************************************************/
		free_dvector(rtr,1,M);
		free_dvector(rti,1,M); 
		free_dvector(p,1,n);
		free_dvector(v,1,n);
//		free(coeff);

		printf("Repeat for another temperature? (Y/N): ");
		scanf("%c", &repeat);

	} while(repeat == 'y' || repeat == 'Y');

	getchar();
}

#undef N
#undef M
#undef b
#undef a
#undef R
