#include <stdio.h>
#include <stdlib.h>

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
    printf("Index of minimum b[i]: %d\n", min_index);
    printf("Minimum b value: %lf\n", b[min_index]);

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
