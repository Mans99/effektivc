#include<stdio.h>
#include<stdlib.h>

struct simplex_t {
    int m;
    int n;
    int* var;
    double ** A;
    double* b;
    double* x;
    double* c;
    double y;

};

int init(struct simplex_t* s,int m,int n,int* var, double** A, double* b, double* x, double* c, double y) {
    int i, k;
    s->b = b;
    s->m = m;
    s->A = A;
    s->c = c;
    s->n = n;
    s->var = var;
    s-> x = x;
    s-> y = y;

    if (s->var == NULL) {
        s->var = calloc(m+m+1, sizeof(int));
        for (i = 0; i < m+n; i=i+1) {
            s->var[i] = i;
        }
    for (k = 0, i = 1; i <m; i++) {
        if (s->b[i] < s->b[k]) {
            k = i;
        }
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
