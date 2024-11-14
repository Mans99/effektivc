#include <stdio.h>
#include <stdlib.h>

void scan_matrix(int m, int n, double** c, double*** a, double** b);
void print_matrix(int m, int n, double* c, double** a, double* b);

int main(int argc, char** argv) {
	int m;
	int n;
	double* c;
	double** a;
	double* b;
	int i;
	int j;


	scanf("%d %d", &m, &n);
	printf("m = %d, n = %d  \n", m, n);
	c = calloc(n, sizeof(double));
	a = calloc(m, sizeof(double*));
	for (i = 0; i < m; i+= 1) {
		a[i] = calloc(n, sizeof(double));
	}
	b = calloc(m, sizeof(double));


	for(i = 0; i < n; i += 1) {
		scanf("%lf", &c[i]);
	}

	for(i = 0; i < m; i += 1) {
		for (j = 0; j < n; j += 1) {
			scanf("%lf", &a[i][j]);
		}
	}

	for(i = 0; i < n; i += 1) {
		scanf("%lf", &b[i]);
	}

	print_matrix(m, n, c, a, b);
	return 0;
}


void print_matrix(int m, int n, double* c, double** a, double* b) {
	int i;
	int j;

	printf("max Z = ");
    for(i = 0; i < n; i +=1){
        if(i > 0){
            printf("%+10.3lf ",c[i]);
        } else {
            printf("%10.3lf", c[i]);
        }
    }
    printf("\n");

    for (i = 0; i < m; i += 1) {
        printf("\t");
        for (j = 0; j < n; j += 1) {
            if (j > 0) {
                printf("%+10.3lf", a[i][j]);
            } else {
                printf("%10.3lf", a[i][j]);
            }
        }
        printf(" \u2264 %8.3lf\n", b[i]);
    }
}
