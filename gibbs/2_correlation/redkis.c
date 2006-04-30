#include "../../nr/gaussj.c"

void redkis (int n, int m, double *x, double *y, FILE *fp, double *yc) {

	double *v, *b, **a, **ai, **A;
	int i, j, k;
			
	printf(     "\n\nNumber of parameters of Redlich-Kister function: %d\n", m);
	fprintf(fp, "\n\nNumber of parameters of Redlich-Kister function: %d\n", m);
	
	v  = dvector(1,2*m-1);
	b  = dvector(1,m);
	a  = dmatrix(1,m,1,m);
	ai = dmatrix(1,m,1,m);
	A  = dmatrix(1,m,1,1);

	for(j=1; j<=(2*m-1); j++) {
		v[j] = 0;
		for (i=1; i<=n; i++)
			v[j] += pow(x[i]*(1 - x[i]), 2)*pow(2*x[i] - 1, j-1);
	}

	for (k=1; k<=m; k++)
		for (i=1, j=k; i<=k, j>=1; i++, j--)
			a[i][j] = v[k];

	for (k=m+1; k<=(2*m-1); k++)
		for (i=(k-m+1), j=m; i<=m, j>=(k-m+1); i++, j--)
			a[i][j] = v[k];

	for (j=1; j<=m; j++) {
		b[j] = 0.0;
		for (i=1; i<=n; i++)
			b[j] += y[i]*x[i]*(1 - x[i])*pow(2*x[i] - 1, j-1);
	}

	for (i=1;i<=m;i++) {
		for (j=1;j<=m;j++)
			ai[i][j] = a[i][j];
		A[i][1] = b[i];
	}

	gaussj (ai, m, A, 1);
	
	/*************************************************************************/
	
	double s, *S, *xj, yj, *Sjy2, Sy2, R2;

	S    = dvector(1,m);
	xj   = dvector(1,m);
	Sjy2 = dvector(1,m);

	// Corrected values
	for (i=1; i<=n; i++) {
		yc[i] = 0.0;
		for(j=1; j<=m; j++)
			yc[i] += x[i]*(1 - x[i])*(A[j][1]*pow(2*x[i] - 1, j-1));
	}
		
	// Standard deviation
	s = 0.0;
	for (i=1; i<=n; i++) 
		s += pow(y[i] - yc[i], 2);
	s = sqrt(s / (double)(n - m) );
	printf(     "\nStandard deviation of the fit: \ns = %3.3lf\n", s);
	fprintf(fp, "\nStandard deviation of the fit: \ns = %3.3lf\n", s);
	
	// Ai uncertainties
	for(i=1; i<=m; i++)
		S[i] = s*sqrt(ai[i][i]);
	
	printf(     "\nFit parameters: A = ai*b, (ai = inverted a):\n");
	fprintf(fp, "\nFit parameters: A = ai*b, (ai = inverted a):\n");
	printf(     "\n%10s %10s\n", "A[i][1]", "S[i]");
	fprintf(fp, "\n%10s %10s\n", "A[i][1]", "S[i]");
	for (i=1; i<=m; i++) {
		printf(     "%10.3lf %10.3lf\n", A[i][1], S[i]);
		fprintf(fp, "%10.3lf %10.3lf\n", A[i][1], S[i]);
	}
	
	// Regresion coefficients
	for (j=1; j<=m; j++) {
		xj[j] = 0.0;
		for (i=1; i<=n; i++)
			xj[j] += pow(x[i], j-1);
		xj[j] = xj[j] / (double)n;
	}
	
	yj = 0.0;
	for (i=1; i<=n; i++)
		yj += y[i];
	yj = yj / (double)n;
	
	for (j=1; j<=m; j++) {
		Sjy2[j] = 0.0;
		for (i=1; i<=n; i++)
			Sjy2[j] += (pow(x[i], j-1) - xj[j])*(y[i] - yj);
		Sjy2[j] = Sjy2[j] / (double)(n-1);
	}
	
	Sy2 = 0.0;
	for (i=1; i<=n; i++)
		Sy2 += pow(y[i] - yj, 2);
	Sy2 = Sy2 / (double)(n-1);
			
	R2 = 0.0;
	for (j=1; j<=m; j++)
		R2 += A[j][1]*Sjy2[j] / Sy2;		
	printf(     "\nRegression coefficient: R2 = %lf\n", R2);
	fprintf(fp, "\nRegression coefficient: R2 = %lf\n", R2);
	
	free_dvector(Sjy2,1,m);
	free_dvector(xj,1,m);
	free_dvector(S,1,m);
	free_dmatrix(A,1,m,1,1);
	free_dmatrix(ai,1,m,1,m);
	free_dmatrix(a,1,m,1,m);
	free_dvector(b,1,m);
	free_dvector(v,1,2*m-1);
}

