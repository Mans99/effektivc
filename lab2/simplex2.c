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

int init(struct simplex_t* s, int m, int n, double** A, double* b, double* c, double* x, double y , int* var);
int select_nonbasic(struct simplex_t* s);
void prepare(struct simplex_t* s,int k);
int initial(struct simplex_t* s,int m,int n, double** A, double* b, double* c, double* x, double y,int* var);
void pivot (struct simplex_t* s, int row, int col);
double xsimplex(int m, int n, double** A, double* b, double* c, double* x, double y, int* var, int h);
int simplex(int m, int n, double ** A, double* b, double* c, double* x, double y);



int init(struct simplex_t* s, int m, int n, double** A, double* b, double* c, double* x, double y , int* var) {
    int i, k;
    s->m = m;
    s->n = n;
    s->A = A;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;
    s->var = NULL;

    // Allocate `var` if it is NULL
    if (s->var == NULL) {
        s->var = calloc(m + n + 1, sizeof(int));
        for (i = 0; i < m + n; i++) {
            s->var[i] = i;
        }
    }

    // Find index `k` of the smallest `b[i]`
    k = 0;
    for (i = 1; i < m; i++) {
        if (b[i] < b[k]) {
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



int initial(struct simplex_t* s,int m,int n, double** A, double* b, double* c, double* x, double y, int* var){
    int i,j,k;
    double w;
    k = init(s, m, n, A, b, c, x, y, var);

     if (b[k] >= 0){
        return 1;
    }
    prepare(s,k);

    n = s->n;

    s->y = xsimplex(m,n,s->A,s->b,s->c,s->x,0,s->var,1);

    for (i = 0; i < m+n; i++){
        if (s->var[i] == m+n-1){
            if (fabs(s->x[i]) > pow(10,-6)){
                free(s->x);
                free(s->c);
                return 0;
            } else {
                break;
            }
        }
        if (i >= n){
            for (j = k = 0; k < n; k++){
                if (fabs(s->A[i-n][k])>fabs(s->A[i-n][j])){
                    j = k;
                }
            }
            pivot(s,i-n,j);
            i=j;
        }
        if (i < n-1){
            k = s->var[i]; 
            s->var[i] = s->var[n-1];
            s->var[n-1]=k;
            for (int k = 0; k < m; k++){
                w = s->A[k][n-1];
                s->A[k][n-1] = s->A[k][i];
                s->A[k][i] = w;
            }
        }
        free(s->c);
        s->c = c;
        s->y = y;
        for (k = n-1; k < n+m-1; k++){
            s->var[k] = s->var[k+1];
        }
        n = s->n = s->n-1;
        double* t = calloc(n, sizeof(double));
        for (k = 0; k < n; k++){
            for (j = 0; j < n; j++){
                if (k == s->var[j]){
                    t[j] = t[j] + s->c[k];
                    break;
                }
            }
            for (j=0; j<m;j++){
                if (s->var[n+j] == k){
                    break;
                }
            }
            s->y = s->y + s->c[k]*s->b[j];
            for (i=0;i<n;i++){
                t[i]=t[i] - s->c[k]*s->A[j][i];
            }
        }
        for (i = 0; i<n;i++){
            s->c[i]=t[i];
        }
        free(t);
        free(s->x);
        return 1;
    }
    return 0;
}



void pivot (struct simplex_t* s, int row, int col) {
    double** a = s->A;
    double* b = s->b;
    double* c = s->c;
    int m = s->m;
    int n = s->n;
    int i,j,t;

    t = s->var[col];
    s->var[col] = s->var[n + row];
    s->var[n+row] = t;
    s->y = s->y + c[col] * b[row]/a[row][col];
    for (i = 0; i < n ; i = i + 1) {
        if (i != col) {
            c[i] = c[i] - c[col] * a[row][i] / a[row][col];
        }
    }
    c[col] = -c[col] / a[row][col];
    for (i = 0; i < m; i = i+1) {
        if (i != row) {
            b[i] = b[i] - a[i][col] * b[row]/ a[row][col];
        }
    }
    for (i = 0; i < m; i = i+1) {
        if(i != row) {
            for (j = 0; j < n; j++) {
                if (j != col) {
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }
    for (i = 0; i < m; i++) {
        if (i != row) {
            a[i][col] = -a[i][col] / a[row][col];
        }
    }
    for (i = 0; i < n; i++) {
        if (i != col) {
            a[row][i] = a[row][i] / a[row][col];
        }
    }
    b[row] = b[row] / a[row][col];
    a[row][col] = 1 / a[row][col];
}



double xsimplex(int m, int n, double** A, double* b, double* c, double* x, double y, int* var, int h) {
    struct simplex_t s;
    int i, row, col;


    if (!initial(&s, m, n, A, b, c, x, y, var)) {
        free(s.var);
        return NAN; 
    }

    while ((col = select_nonbasic(&s)) >= 0) {
        row = -1;

        for (i = 0; i < m; i++) {
            if (A[i][col] > pow(10,-6)) {
                if (row < 0 || (b[i] / A[i][col] < b[row] / A[row][col])) {
                    row = i;
                }
            }
        }

        if (row < 0) {
            free(s.var);
            return INFINITY; 
        }

        pivot(&s, row, col);
    }


    if (h == 0) {
        for (i = 0; i < n; i++) {
            if (s.var[i] < n) {
                x[s->var[i]] = 0;
            }
        }
        for (i = 0; i < m; i++) {
            if (s->var[n + i] < n) {
                x[s->var[n + i]] = s->b[i];
            }
        }
    } else {
        for (i = 0; i < n; i++) {
            x[i] = 0;
        }
        for (i = n; i < n + m; i++) {
            x[i] = s->b[i - n];
        }
    }

    double result = s.y;
    free(s.var);
    return result;
}


int simplex(int m, int n, double ** A, double* b, double* c, double* x, double y) {
    return xsimplex(m,n,A,b,c,x,y, NULL,0);
}




int main(int argc, char** argv) {
    int m, n;
    double* c;
    double** a;
    double* b;
    int i, j;
    int* var; 
    double y = 0;
    double* x = calloc((n+1),sizeof(double)); 

    scanf("%d %d", &m, &n);
    printf("m = %d, n = %d\n", m, n);

    c = calloc(n, sizeof(double));
    for (i = 0; i < n; i++) {
        printf("%lf", c[i]);
    }
    
    a = calloc(m, sizeof(double*));
    for (i = 0; i < m; i++) {
        a[i] = calloc(n + 1, sizeof(double));
    }
    b = calloc(m, sizeof(double));

    for (i = 0; i < n; i++) {
        scanf("%lf", &c[i]);
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            scanf("%lf", &a[i][j]);
        }
    }

    for (i = 0; i < m; i++) {
        scanf("%lf", &b[i]);
    }

    int min_index = simplex(m, n, a, b, c, x, y);

    printf("Index of minimum b[i]: %d\n", min_index);

    free(c);
    for (i = 0; i < m; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);

    return 0;
}




//step: 1 rad
//next: nästa metod
//up	Move up the stack trace to the calling function.
//down	Move down the stack trace to the called function.
//display	show value of variable 