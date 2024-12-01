#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <set>


struct node_t
{
    int m; /* Constraints. */
    int n; /* Decision variables. */
    int k; /* Parent branches on xk. */
    int h; /* Branch on xh. */
    double xh; /* xh. */
    double ak; /* Parent ak. */
    double bk; /* Parent bk. */
    double* min; /* Lower bounds. */
    double* max; /* Upper bounds. */
    double** a; /* A. */
    double* b; /* b. */
    double* x; /* x. */
    double* c; /* c. */
    double z; /* z. */
};



void bound(struct node_t* p,int h, double zp, double* x);
int branch(struct node_t* q, double z);
void succ(struct node_t* p,int h, int m, int n, double** a, double* b,double* c,int k,double ak,double bk, double* zp, double* x);
int intopt(int m, int n, double** a, double* b, double* c, double* x);

void bound(struct node_t* p,int h, double* zp, double* x) {
    if (p.z > *zp) {
        *zp = p.z;
        memcpy(x, p.x, sizeof(double) * p.n);
        

    } 
}

int branch(struct node_t* q, double z) {
    double min, max;
    if (q.z < z) {
        return 0;
    }
    for(h=0; h < q.n; h = h+1) {
        if (!is_integer(&q.x[h])) {
            if(q.min[h] = -INFINITY) {
                min = 0;
            } else {
                min = q.min[h];
            }
            max = q.max[h];
            if(floor(q.x[h]) < min || ceil(q.x[h] > max)) {
                continue;
            }
            q.h = h;
            q.xh = q.x[h];

            for (int i = 0; i < m; i++) {
                free(q.a[i]);
            }
            free(q.a);
            free(q.b);
            free(q.c);
            free(q.x);
            return 1;
        }
    return 0;
    }
} 


void succ(struct node_t* p,int h, int m, int n, double** a, double* b,double* c,int k,double ak,double bk, double* zp, double* x) {
    struct node_t* q = extend(p,m,n,a,b,c,k,ak,bk);
    if(q = null) {
        return;
    }
    q.z = simplex(q.m, q.n, q.a, q.b, q.c, q.x, 0);
    if isfinite(q.z) {
        if (integer(q)) {
            bound(q,h,zp,x);
        } else if (branch(q, *zp)) {
            add(h, q);
            return;
        }
    }
    free(q);
}

int intopt(int m, int n, double** a, double* b, double* c, double* x) {
    struct node_t p = initial_node(m, n, a, b, c);
    set h = {p};
    double z = -INFINITY;
    p.z = simplex(p.m, p.n, p.a, p.b, p.c, p.x, 0);
    if (integer(p) || !isfinite(p.z)) {
        z = p.z;
        if (integer(p)) {
            x = p.x;
        }
        free(p);
        return z;
    }
    branch(p,z);
    while (h != NULL) {
        p = h;
        succ(p, h, m, n, a, b, c, p.h, 1, floor(p.xh), &z, x);
        succ(p, h, m, n, a, b, c, p.h, -1, -ceil(p.xh), &z, x);
        free(p);
    }
    if(z = -INFINITY) {
        return NAN;
    } else {
        return z;
    }
}


