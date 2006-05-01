void vecfunc(int l, double f[], double fvec[]) {
    fvec[1] = f[1] - ((x1*p1o) / f[2]) * exp((m1 - ((f[2] - p1o)*(B11 - v1o) + d12*pow(1 - f[1], 2)*f[2])*0.001) / R*T);
	fvec[2] = (1 - f[1]) - (((1 - x1)*p2o)/f[2]) * exp((m2 - ((f[2] - p2o)*(B22 - v2o) + d12*f[1]*f[1]*f[2])*0.001) / R*T);
}
