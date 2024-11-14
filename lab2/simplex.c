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

    if (*s->var == NULL) {
        *s->var = calloc(m+m+1, sizeof(int));
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

