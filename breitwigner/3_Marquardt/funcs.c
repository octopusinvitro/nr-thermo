#include <math.h>

void funcs(double x, double a[], double *y, double dyda[], int N) {

	double A, B, den;

	A	= pow(x - a[2], 2);
	B	= pow(a[3], 2) / 4.0;
	den = A + B;
	N	= 3;

	*y = a[1]/den;

	dyda[1] = 1 / den;
	dyda[2] = 2.0*a[1]*(x - a[2]) / pow(den, 2);
	dyda[3] = - (a[1]*a[3]) / (2.0*pow(den, 2));
}

