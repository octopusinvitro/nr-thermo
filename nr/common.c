#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/***************************************************
 * COUNT THE NUMBER OF LINES IN THE DATA FILE 
 * (Needed for NR vector/matrix declarations)
 ***************************************************/
int numOfLines(char *filename) {
	FILE *fp;
	int n = 0;
	char line[100];
	
	if ( (fp = fopen(filename, "r")) == NULL ) {
		printf("\nError opening file\n");
		exit(1);
	}
	
	do {
		fgets(line, 100, fp);
		if ( strlen(line) > 2 )
			n++;
	} while(!feof(fp));
	
	fclose(fp);
	return n;
}

/***************************************************
 * READ DATA FILE
 * (Format: one column of data)
 ***************************************************/
void read1col(int n, double *x, char *filename) {
	FILE *fp;

	if ( (fp=fopen(filename, "r")) == NULL ) {
		printf("\nError opening file\n");
		exit(1);
	}

	else {
		int i;
		for (i=1; i<=n; i++) {
			fscanf(fp,"%lf",&x[i]);
		}
	}

	fclose(fp);
}

/***************************************************
 * READ DATA FILE
 * (Format: two columns of data separated by a space)
 ***************************************************/
void read2col(int n, double *x, double *y, char *filename) {
	FILE *fp;

	if ( (fp=fopen(filename, "r")) == NULL ) {
		printf("\nError opening file\n");
		exit(1);
	}

	else {
		int i;
		for (i=1; i<=n; i++) {
			fscanf(fp, "%lf %lf", &x[i], &y[i]);
		}
	}
	
	fclose(fp);
}

/***************************************************
 * READ DATA FILE
 * (Format: three columns of data separated by a space)
 ***************************************************/
void read3col(int n, double *x, double *y, double *sig, char *filename) {
	FILE *fp;

	if ( (fp=fopen(filename, "r")) == NULL ) {
		printf("\nError opening file\n");
		exit(1);
	}

	else {
		int i;
		for (i=1; i<=n; i++) {
			fscanf(fp, "%lf %lf %lf", &x[i], &y[i], &sig[i]);
		}
	}
	
	fclose(fp);
}

/***************************************************
 * GET THE MAXIMUM VALUE IN A VECTOR OF DOUBLES 
 * (Since Math lacks a method for that)
 ***************************************************/
int maximum (double y[], int n) {
	int i, m;
	double big = -1000000000;

	for (i=1; i<=n; i++) {
		if (y[i] > big) {
			big = y[i];
			m   = i;
		}
	}

	return m;
}
