#include <stdio.h>
#include <stdlib.h>

void scan_matrix(int m, int n, double** c, double*** a, double** b);
void print_matrix(int m, int n, double** c, double*** a, double** b);
double** make_matrix(int m, int n);

int main(int argc, char** argv) {
	int m;
	int n;
	double* c;
	double** a;
	double* b;


	scanf("%d %d", &m, &n);
	printf("m = %d, n = %d  \n", m, n);
	c = calloc(n, sizeof(double));
	a = make_matrix(m, n);
	b = calloc(m, sizeof(double));


	scan_matrix(m, n, &c, &a, &b);
	print_matrix(m, n, &c, &a, &b);
	return 0;
}


void scan_matrix(int m, int n, double** c, double*** a, double** b) {
	int i;
	int j;
	double ** s = *a;

	
	
	for(i = 0; i < n; i += 1) {
		scanf("%lf", c[i]);
	}

	for(i = 0; i < m; i += 1) {
		for (j = 0; j < n; j += 1) {
			scanf("%lf", &s[i][j]);
		}
	}

	for(i = 0; i < n; i += 1) {
		scanf("%lf", b[i]);
	}

}

void print_matrix(int m, int n, double** c, double*** a, double** b) {
	int i;
	int j;
	double ** s = *a;
	printf("C = \n");
	

	for(i = 0; i < n; i += 1) {
		printf("%10.3lf\n", *c[i]);
	}

	printf("A = \n");

	for(i = 0; i < m; i += 1) {
		for (j = 0; j < n; j += 1) {
			printf("%+10.3lf\n", s[i][j]);
		}
	}

	printf("B = \n");

	for(i = 0; i < m - 1; i += 1) {
		printf("%10.3lf ", *b[i]);
	}

	printf("%10.3lf \n", *b[n -1]);

}


double** make_matrix(int m, int n) {
	int i;
	double ** a;
	a = calloc(m, sizeof(double*));
	for (i = 0; i < m; i+= 1) {
		a[i] = calloc(n, sizeof(double));
	}
	return a;
}
