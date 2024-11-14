#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct simplex_t {
    int m;
    int n;
    int* var;
    double **A;
    double* b;
    double* x;
    double* c;
    double y;
};

int init(struct simplex_t* s, int m, int n, int* var, double** A, double* b, double* x, double* c, double y) {
    int i, k;
    s->m = m;
    s->n = n;
    s->A = A;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;

    // Allocate `var` if it is NULL
    if (s->var == NULL) {
        s->var = calloc(m + n, sizeof(int));
        for (i = 0; i < m + n; i++) {
            s->var[i] = i;
        }
    }

    // Find index `k` of the smallest `b[i]`
    k = 0;
    for (i = 1; i < m; i++) {
        if (s->b[i] < s->b[k]) {
            k = i;
        }
    }
    return k;
}


int select_nonbasic(struct simplex_t* s){
    int i;
    for (int i = 0; i < s->n; i++){
        if (s->c[i] > pow(10, -6)){
            return i;
        }
    }
    return -1;
}

/*

void prepare(struct simplex_t* s,int k){
    int m = s->m;
    int n = s->n;
    int i;

    for (int i = m + n; i > n; i--){
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m+n;
    n = n + 1;
    for (int i=0; i < m; i++){
        s->A[i][n-1] = s->A[i][n-1] - 1;
    }
    s->x = calloc(m+n, sizeof(double));
    s->c = calloc(n, sizeof(double));
    s->c[n-1] = -1;
    s->n=n;
    pivot(s,k,n-1);
}

int initial(struct simplex_t* s,int m,int n,int* var, double** A, double* b, double* x, double* c, double y){
    int i,j,k;
    double w;
    k = init(&s, m, n, var, A, b, x, c, y);

    if (b[k] >= 0){
        return 1;
    }
    prepare(&s,k);

    n = s->n;
    s->y = xsimplex(m,n,s->A,s->b,s->c,s->x,0,s->var,1);

    for (i = 0; i < m+n+1; i++){
        if (s->var[i] == m+n+1){
            if (abs(s->x[i]) > pow(10,-6)){
                free(s->x);
                free(s->c);
                return 0;
            } else {
                break;
            }
        }
        if (i >= n){
            for (j = k = 0; k < n; k++){
                if (abs(s->A[i-n][k])>abs(s->A[i-n][j])){
                    j = k;
                }
            }
            pivot(s,i-n,j);
            i=j;
        }

    }


}

*/


int main(int argc, char** argv) {
    int m, n;
    double* c;
    double** a;
    double* b;
    int i, j;
    struct simplex_t simplex;
    int* var = NULL;  // Initialize as NULL to be allocated in `init()`
    double y = 0;
    double* x = NULL;  // This can remain NULL if not used in the current logic

    // Read m and n
    scanf("%d %d", &m, &n);
    printf("m = %d, n = %d\n", m, n);

    // Allocate memory for c, a, and b
    c = calloc(n, sizeof(double));
    a = calloc(m, sizeof(double*));
    for (i = 0; i < m; i++) {
        a[i] = calloc(n, sizeof(double));
    }
    b = calloc(m, sizeof(double));

    // Read coefficients of the objective function c
    for (i = 0; i < n; i++) {
        scanf("%lf", &c[i]);
    }

    // Read matrix A
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            scanf("%lf", &a[i][j]);
        }
    }

    // Read vector b
    for (i = 0; i < m; i++) {
        scanf("%lf", &b[i]);
    }

    // Initialize the simplex structure
    int min_index = init(&simplex, m, n, var, a, b, x, c, y);
    int nonbasic = select_nonbasic(&simplex);

    printf("Index of minimum b[i]: %d\n", min_index);
    printf("Minimum b value: %lf\n", b[min_index]);
    printf("Nonbasic value: %lf", c[nonbasic]);

    // Free allocated memory
    free(c);
    for (i = 0; i < m; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(simplex.var);

    return 0;
}