#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*struct note_t
initial_node
extend
is_integer
integer*/

struct node_t{
    int m;
    int n;
    int k;
    int h;
    double xh;
    double ak;
    double bk;
    double* min;
    double* max;
    double** a;
    double* b;
    double* x;
    double* c;
    double z;
};

struct node_t initial_node(int m, int n, double** a, double* b, double* c){
    int i, j;
    struct node_t p;
    p.a = calloc(m+1, sizeof(double*));
    for (i = 0; i < m+1; i++){
        p.a[i] = calloc(n+1, sizeof(double));
    }
    p.b = calloc(m+1, sizeof(double));
    p.c = calloc(n+1, sizeof(double));
    p.x = calloc(n+1, sizeof(double));
    p.min = calloc(n, sizeof(double));
    p.max = calloc(n, sizeof(double));
    p.m = m;
    p.n = n;

    //copy a, b, and c parameters to p ??

    for (i = 0; i < m+1; i++){
        p.b[i] = b[i];
        for (j = 0; j < n+1; j++){
            p.a[i][j] = a[i][j];
 
        }
    }

    for (i = 0; i < n+1; i++) {
        p.c[i] = c[i];
    }


    for (i = 0; i < n; i++){
        p.min[i] = -INFINITY;
        p.max[i] = INFINITY;
    }
    return p;
}

struct node_t extend(struct node_t p,int m, int n, double** a, double* b, double* c, int k, double ak, double bk){
    struct node_t q;
    int i,j;
    q.k = k;
    q.ak = ak;
    q.bk = bk;
    if (ak > 0 && p.max[k] < INFINITY){
        q.m = p.m;
    } else if (ak < 0 && p.min[k]>0){
        q.m = p.m;
    } else {
        q.m = p.m + 1;
    }
    q.n = p.n;
    q.h = -1;
    q.a = malloc((q.m+1) * sizeof(double*));
    for (i = 0; i < q.m+1; i++){
        q.a[i] = malloc((q.n+1) *sizeof(double));
    }
    q.b = malloc((q.m+1)*sizeof(double));
    q.c = malloc((q.n+1)*sizeof(double));
    q.x = malloc((q.n+1)*sizeof(double));
    q.min = malloc((n)*sizeof(double));
    q.max = malloc((n)*sizeof(double));

    for (i = 0; i < n; i++){
        q.min[i] = p.min[i];
        q.max[i] = p.max[i];   
    }

    for (i = 0; i < m; i++){
        q.b[i] = b[i];
        for (j = 0; j < n+1; j++){
            q.a[i][j] = a[i][j];
        }
    }
    for (i = 0; i < n+1; i++){
        q.c[i] = c[i];
    }

    if (ak > 0){
        if (q.max[k] == INFINITY || bk < q.max[k]){
            q.max[k] = bk;
        }
    } else if (q.min[k] == -INFINITY || -bk > q.min[k]){
            q.min[k] = -bk;
    }
    for (i = m, j = 0; j < n; j++){
        if (q.min[j] > -INFINITY){
            q.a[i][j] = -1;
            q.b[i] = -q.min[j];
            i++;
        } if (q.max[j] < INFINITY) {
            q.a[i][j] = 1;
            q.b[i] = q.max[j];
            i++;
        }
    }
    return q;
    }

    int is_integer(double* xp){
        double x = *xp;
        double r = lround(x);
        if (abs(r-x) < pow(10,-6)){
            *xp = r;
            return 1;
        } else {
            return 0;
        }
    }

    int integer(struct node_t p){
        int i;
        for (i = 0; i < p.n; i++){
            if (!is_integer(&p.x[i])){
                return 0;
            }
        }
        return 1;
    }

