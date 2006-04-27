/**********************************************************************************
 * MARQUARDT: Program that obtains the parameters of a (non-linear) user function, 
 *            using the method of least-squares and the method of
 *            Levemberg-Marquardt.
 *
 * APPLICATION: To obtain the cross-section according to Breit-Wigner formula.
 **********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "funcs.c"
#include "../../nr/nrutil.c"
#include "../../nr/covsrt.c"
#include "../../nr/gaussj.c"
#include "../../nr/mrqcof.c"
#include "../../nr/mrqmin.c"
#include "../../nr/common.c"

int main(void) {

	int numOfLines(char *filename), maximum (double y[], int N);
	void read3col(int n, double *x, double *y, double *sig, char *filename);
	void funcs(double x, double a[], double *y, double dyda[], int N);
	void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[], 
                int ma, double **covar, double **alpha, double *chisq, 
                void (*funcs)(double, double [], double *, double [], int), double *alamda);

	int M, N, *ia, i, j, k, iter = 0, n = 0;
	double *x, *y, *sig, sum = 0.0, *a, **covar, **alpha, chisq, alamda, ochisq, *ymod, *dyda;

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
	M = numOfLines(name);
	printf(     "Number of points = %d\n", M);
	fprintf(fp, "Number of points = %d\n", M); 
	getchar();

	printf(     "Number of parameters of the function: ");
	scanf("%d", &N);
	fprintf(fp, "Number of parameters of the function = %d\n\n", N); 
	getchar();
	 
	/***************************************************
	 * NR vector and matrix declaration
	 ***************************************************/
	x     = dvector(1, M);
	y     = dvector(1, M);
	sig   = dvector(1, M);
	a     = dvector(1, N); 
	ymod  = dvector(1, M);
	dyda  = dvector(1, N);
	ia    = ivector(1, N);
	covar = dmatrix(1, N, 1, N);
	alpha = dmatrix(1, N, 1, N);
	 
	/***************************************************
	 * Filling in vectors and matrices
	 ***************************************************/
	read3col(M, x, y, sig, name);
	printf( "\n\nVectors x and y:\n");
	fprintf(fp, "Vectors x and y:\n");
	for(i=1; i<=M; i++) {
		printf(     "\n E[%d] = %3.3lf \t Sig[%d] = %3.3lf \t dSigma[%d] = %3.3lf", i, x[i], i, y[i], i, sig[i]);
		fprintf(fp, "\n E[%d] = %3.3lf \t Sig[%d] = %3.3lf \t dSigma[%d] = %3.3lf", i, x[i], i, y[i], i, sig[i]);
	}
	printf(     "\n\n\n");
	fprintf(fp, "\n\n\n");

	/***************************************************
	 * Initial values for the parameter vectors
	 * For example:
	 *    a[1] = 70000.0;
	 *    a[2] = 75.0;
	 *    a[3] = 60.0;
	 ***************************************************/
	for (i=1; i<=N; i++) {
		printf("Enter initial value of parameter %d: ", i);
		scanf("%lf", &a[i]);
	}

	// Fit all
	for (i=1; i<=N; i++)
		ia[i] = 1;

	chisq  = 0.0;
	alamda = -1;
	
	printf(     "\nParameters (initial values):\n");
	fprintf(fp, "\nParameters (initial values):\n");
	for (i=1; i<=N; i++) {
		printf(     " a[%d] = %3.3lf \n", i, a[i]);
		fprintf(fp, " a[%d] = %3.3lf \n", i, a[i]);
	}

	getchar();
		
	/***************************************************
	 * ITERATION LOOP
	 ***************************************************/
	char answer;
	do {
		printf("\nHow many iterations?: ");
		scanf("%d", &k);

		for(i=1; i<=k; i++) {
			iter++;
			ochisq = chisq;

			mrqmin(x, y, sig, M, a, ia, N, covar, alpha, &chisq, funcs, &alamda);

			if (chisq > ochisq)
				n=0;
			else if (fabs(ochisq - chisq) < 0.1)
				n++;

			if(n >= 4)
				break;
		}

		if(n < 4) {
			printf("\nIterate again? (Y/N)?: ");
			scanf("%c", &answer);
		}
	} while (answer == 's' || answer == 'S');
	
	// Last iterations
	alamda = 0.0;
	mrqmin(x, y, sig, M, a, ia, N, covar, alpha, &chisq, funcs, &alamda);
	iter++;
	 
	/***************************************************
	 * PRINTING RESULTS
	 ***************************************************/
	printf(     "\nExit on iter = %d \n", iter);
	fprintf(fp, "\nExit on iter = %d \n", iter);
	printf(     "\nParameters (after calculations):\n");
	fprintf(fp, "\nParameters (after calculations):\n");
	for (i=1; i<=N; i++) {
		printf(     " a[%d] = %3.3lf \t da[%d]: %3.3lf \n", i, a[i], i, sqrt(covar[i][i]));
		fprintf(fp, " a[%d] = %3.3lf \t da[%d]: %3.3lf \n", i, a[i], i, sqrt(covar[i][i]));
	}
	printf(     "\nChi squared = %3.3lf \n", chisq);
	fprintf(fp, "\nChi squared = %3.3lf \n", chisq);
	printf(     "n = %d\n\n", n);
	fprintf(fp, "n = %d\n\n", n);
	getchar();

	printf(     "\nCovariance matrix:\n\n");
	fprintf(fp, "\nCovariance matrix:\n\n");
	for (i=1; i<=N; i++) {
		for (j=1; j<=N; j++) {
			printf(     " covar[%d][%d] = %3.3lf \t", i, j, covar[i][j]);
			fprintf(fp, " covar[%d][%d] = %3.3lf \t", i, j, covar[i][j]);
		}
		printf(     "\n");
		fprintf(fp, "\n");
	}
	printf(     "\n");
	fprintf(fp, "\n");

	printf(     "\nAlpha matrix:\n\n");
	fprintf(fp, "\nAlpha matrix:\n\n");
	for (i=1; i<=N; i++) {
		for (j=1; j<=N; j++) {
			printf(     " alpha[%d][%d] = %3.3lf \t", i, j, alpha[i][j]);
			fprintf(fp, " alpha[%d][%d] = %3.3lf \t", i, j, alpha[i][j]);
		}
		printf(     "\n");
		fprintf(fp, "\n");
	}

	printf(     "\nData points recalculated from the fit:\n");
	fprintf(fp, "\nData points recalculated from the fit:\n");
	for (i=1; i<=M; i++) {

		(*funcs)(x[i], a, &ymod[i], dyda, N);
	
		sum = 0.0;
		for (j=1; j<=N; j++)
			sum += pow(dyda[j]*sqrt(covar[j][j]), 2);

		printf(     "\nSigma[%d] = %3.3lf\t\tdSigma[%d] = %3.3lf", i, ymod[i], i, sqrt(sum));
		fprintf(fp, "\nSigma[%d] = %3.3lf\t\tdSigma[%d] = %3.3lf", i, ymod[i], i, sqrt(sum));
	}
		 
	/***************************************************
	 * Free memory
	 ***************************************************/
	free_dmatrix(alpha, 1, N, 1, N);
	free_dmatrix(covar, 1, N, 1, N);
	free_ivector(ia,    1, N);
	free_dvector(a,     1, N);
	free_dvector(dyda,  1, N);
	free_dvector(ymod,  1, M);
	free_dvector(sig,   1, M);
	free_dvector(y,     1, M);
	free_dvector(x,     1, N);

	printf("\n\n***************************************************");
	printf("\n\n               END OF CALCULATIONS");
	printf("\n\n***************************************************\n\n");
	return 0;
}
