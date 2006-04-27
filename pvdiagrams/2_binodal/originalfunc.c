double func(double p, double T) {

	int i, j;
	double *rti, *rtr, aux, r, s, t, area1, area2, c[N];	
	bool realroots = true;
	void zrhqr(double c[], int m, double rtr[], double rti[]);

	c[0] = p;
	c[1] = -(p*b+R*T);
	c[2] = a;
	c[3] = -a*b;
	
	rti = dvector(1,M);
	rtr = dvector(1,M);	
	
	zrhqr(c, M, rtr, rti);
	
	if( (rti[1] == 0.0) && 
		(rti[2] == 0.0) && 
		(rti[3] == 0.0)) {
		for (i=1; i<M; i++) {
			for (j=i+1; j<=M; j++) {
				if(rtr[i]>rtr[j]) {
					aux=rtr[i];
					rtr[i]=rtr[j];
					rtr[i]=aux;
				} 
			}
		}

		r = rtr[1];
		s = rtr[2];
		t = rtr[3];

		area1 = p*(s-r) - 
				R*T*log(fabs((s - b) / (r - b))) - 
				a*((1 / s) - (1 / r));

		area2 = R*T*log(fabs((t - b) / (s - b))) + 
				a*((1 / t) - (1 / s)) - 
				p*(t - s);

//printf("\n area1: %lf, area2: %lf, area1+area2: %lf\n",area1,area2,area1+area2);
	}

	else {
		printf("\n Complex roots \n");
		realroots = false;
	}

	free_dvector(rtr,1,M);
	free_dvector(rti,1,M);
	free(c);

	if (realroots)
		return area1 + area2;
	return 0.0;
}
