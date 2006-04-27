double func(double p, double T) {

	int i, j;
	double *rti, *rtr, aux, area, c[N];
	bool realroots = true;
	void zrhqr(double c[], int m, double rtr[], double rti[]);

	c[0] = -a*b;
	c[1] = a;
	c[2] = -(p*b+R*T);
	c[3] = p;
	
	rti = dvector(1,M);
	rtr = dvector(1,M);	
	
	zrhqr(c, M, rtr, rti);
	
	/* We're only interested in real roots */
	if( (rti[1] == 0.0) && 
		(rti[2] == 0.0) && 
		(rti[3] == 0.0) ) {
		for (i=1; i<M; i++) {
			for (j=i+1; j<=M; j++) {
				if(rtr[i] > rtr[j]) {
					aux    = rtr[i];
					rtr[i] = rtr[j];
					rtr[i] = aux;
				} 
			} 
		}

		area =  R*T*log(fabs((rtr[3] - b) / (rtr[1] - b))) + 
				a*((1 / rtr[3]) - (1 / rtr[1])) - 
				p*(rtr[3] - rtr[1]);
	}
		 
	else {
		printf("\n Complex roots");
		realroots = false;
	}

	free_dvector(rtr,1,M);
	free_dvector(rti,1,M);
//	free(c);

	if (realroots == true)
		return area;
	else {
		return 0.0;
	}
}
