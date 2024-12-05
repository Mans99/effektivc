#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

struct Node {
    struct node_t* p;
    struct Node* next;
};

int init(struct simplex_t* s, int m, int n, double** A, double* b, double* c, double* x, double y , int* var);
int select_nonbasic(struct simplex_t* s);
void prepare(struct simplex_t* s,int k);
int initial(struct simplex_t* s,int m,int n, double** A, double* b, double* c, double* x, double y,int* var);
void pivot (struct simplex_t* s, int row, int col);
double xsimplex(int m, int n, double** A, double* b, double* c, double* x, double y, int* var, int h);
double simplex(int m, int n, double ** A, double* b, double* c, double* x, double y);
double intopt(int m, int n, double** a, double* b, double* c, double* x);
void bound(struct node_t* p,struct Node** h, double* zp, double* x);
int branch(struct node_t* q, double z);
void succ(struct node_t* p,struct Node**  h, int m, int n, double** a, double* b,double* c,int k,double ak,double bk, double* zp, double* x);

int glob = 0;


int init(struct simplex_t* s, int m, int n, double** A, double* b, double* c, double* x, double y , int* var) {
    int i, k;
    s->m = m;
    s->n = n;
    s->A = A;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;
    s->var = var;

    
    

    // Allocate `var` if it is NULL
    if (s->var == NULL) {
        s->var = malloc((m + n + 1) * sizeof(int));
        for (i = 0; i < m + n; i++) {
            s->var[i] = i;
        }
        

    // Find index `k` of the smallest `b[i]`
    k = 0;
    for (i = 1; i < m; i++) {
        if (b[i] < b[k]) {
            k = i;
        }
    }
    return k;
    } else {
    return 0;
    }
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
    glob++;
    int m = s->m;
    int n = s->n;
    int i;


    for (int i = m + n; i > n; i--){
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m+n;
    n = n + 1;
    for (int i=0; i < m; i++){
        s->A[i][n-1] = -1;
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
        } else {

        }
        free(s->c);
        s->c = c;
        s->y = y;
        for (k = n-1; k < n+m-1; k++){
            s->var[k] = s->var[k+1];
        }
        n = s->n - 1;
        s->n = s->n - 1;
        double* t = calloc(n, sizeof(double));
        for (k = 0; k < n; k++){
            for (j = 0; j < n; j++){
                if (k == s->var[j]){
                    t[j] = t[j] + s->c[k];
                    goto next_k;
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
            next_k:;
        }
        for (i = 0; i<n;i++){
            s->c[i]=t[i];
        }
        free(t);
        free(s->x);
        return 1;
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
                x[s.var[i]] = 0;
            }
        }

        for (i = 0; i < m; i++) {
            if (s.var[n + i] < n) {
                x[s.var[n + i]] = s.b[i];
            }
        }
        free(s.var);
    } else {
        for (i = 0; i < n; i++) {
            x[i] = 0;
        }
        for (i = n; i < n + m; i++) {
            x[i] = b[i - n];
        }
    }

    double result = s.y;
    
    return result;
}


double simplex(int m, int n, double ** A, double* b, double* c, double* x, double y) {
    return xsimplex(m,n,A,b,c,x,y, NULL,0);
}


void free_node(struct node_t *node) {
  free(node->min);
  free(node->max);
  free(node->b);
  free(node->x);
  free(node->c);
  for (int i = 0; i < node->m + 1; i++) {
    free(node->a[i]);
  }
  free(node->a);
  free(node);
}



struct Node* create_node(struct node_t* p) {
    struct Node* n = (struct Node*)malloc(sizeof(struct Node));
    n->next = NULL;
    n->p = p;
    return n;
}

void push(struct Node** head, struct node_t* p) {
    struct Node* n = create_node(p);
    n->next = *head;
    *head = n;
}

struct node_t* pop(struct Node** head) {
    if (*head == NULL) {
        return NULL;
    }
    struct Node* n = *head;
    *head = n->next;
    struct node_t* p = n->p;
    free(n);
    return p;
} 


void free_stack(struct Node** head) {
    struct Node* current = *head;
    while (current != NULL) {
        struct Node* temp = current;
        current = current->next;
        free_node(temp->p);  // Free the current node
    }
}

struct node_t* initial_node(int m, int n, double** a, double* b, double* c){
    int i, j;
    struct node_t* p = (struct node_t*)calloc(1, sizeof(struct node_t));
    p->a = calloc(m+1, sizeof(double*));
    for (i = 0; i < m+1; i++){
        p->a[i] = calloc(n+1, sizeof(double));
    }
    p->b = calloc(m+1, sizeof(double));
    p->c = calloc(n+1, sizeof(double));
    p->x = calloc(n+1, sizeof(double));
    p->min = calloc(n, sizeof(double));
    p->max = calloc(n, sizeof(double));
    p->m = m;
    p->n = n;

    for (i = 0; i < m; i++){
        memcpy(p->a[i], a[i], n * sizeof(double));
    }

    memcpy(p->b, b, m * sizeof(double));
    memcpy(p->c, c, n * sizeof(double));


    for (i = 0; i < n; i++){
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY;
    }
    return p;
}

struct node_t* extend(struct node_t* p,int m, int n, double** a, double* b, double* c, int k, double ak, double bk){
    struct node_t* q = (struct node_t*)malloc(sizeof(struct node_t));
    int i,j;
    q->k = k;
    q->ak = ak;
    q->bk = bk;
    if (ak > 0 && p->max[k] < INFINITY){
        q->m = p->m;
    } else if (ak < 0 && p->min[k]>0){
        q->m = p->m;
    } else {
        q->m = p->m + 1;
    }
    q->n = p->n;
    q->h = -1;
    q->a = malloc((q->m+1) * sizeof(double*));
    for (i = 0; i < q->m+1; i++){
        q->a[i] = calloc((q->n+1),sizeof(double));
    }
    q->b = malloc((q->m+1)*sizeof(double));
    q->c = malloc((q->n+1)*sizeof(double));
    q->x = malloc((q->n+1)*sizeof(double));
    q->min = malloc((n)*sizeof(double));
    q->max = malloc((n)*sizeof(double));

    memcpy(q->min, p->min, n * sizeof(double));
    memcpy(q->max, p->max, n * sizeof(double));
    for (i = 0; i < m; i++) {
        memcpy(q->a[i], a[i], n * sizeof(double));
    }
    memcpy(q->b, b, m * sizeof(double));
    memcpy(q->c, c, n * sizeof(double));

    if (ak > 0){
        if (q->max[k] == INFINITY || bk < q->max[k]){
            q->max[k] = bk;
        }
    } else if (q->min[k] == -INFINITY || -bk > q->min[k]){
            q->min[k] = -bk;
    }
    for (i = m, j = 0; j < n; j++){
        if (q->min[j] > -INFINITY){
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i++;
        } if (q->max[j] < INFINITY) {
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i++;
        }
    }
    return q;
    }

    int is_integer(double* xp){
        double x = *xp;
        double r = lround(x);
        if (fabs(r-x) < pow(10,-6)){
            *xp = r;
            return 1;
        } else {
            return 0;
        }
    }

    int integer(struct node_t* p){
        int i;
        for (i = 0; i < p->n; i++){
            if (!is_integer(&(p->x[i]))){
                return 0;
            }
        }
        return 1;
    }




void bound(struct node_t* p,struct Node** h, double* zp, double* x) {
    if (p->z > *zp) {
        *zp = p->z;
        memcpy(x, p->x, sizeof(double) * p->n + 1);
        //*zp = p.z
        //copy each element of p.x to x // save best x
        //remove and delete all nodes q in h with q.z < p.z

        struct Node* current = *h;
        struct Node* prev = NULL;

        while(current != NULL) {
            if (current -> p -> z < p->z) {
                if(prev != NULL) {
                    prev->next = current->next;
                    free_node(current->p);
                    current->next = NULL;
                    free(current);
                    current = prev->next;
                } else {
                    prev = current;
                    current = current -> next;
                    free_node(prev->p);
                    prev->next = NULL;
                    free(prev);
                    *h = current;
                }
            } else {
                prev = current;
                current = current->next;
            }
        }
    }
}

int branch(struct node_t* q, double z) {
    double min, max;
    if (q->z < z) {
        return 0;
    }
    for(int h=0; h < q->n; h = h+1) {
        if (!is_integer(&(q->x[h]))) {
            if(q->min[h] == -INFINITY) {
                min = 0;
            } else {
                min = q->min[h];
            }
            max = q->max[h];
            if(floor(q->x[h]) < min || ceil(q->x[h] > max)) {
                continue;
            }
            q->h = h;
            q->xh = q->x[h];

            // for (int i = 0; i < q->m; i++) {
            //     free(q->a[i]);
            // }
            // free(q->a);
            // free(q->b);
            // free(q->c);
            // free(q->x);
    
            return 1;
        }
    }
    return 0;
} 



void succ(struct node_t* p,struct Node** h, int m, int n, double** a, double* b,double* c,int k,double ak,double bk, double* zp, double* x) {
    struct node_t* q = extend(p,m,n,a,b,c,k,ak,bk);
    if(q == NULL) {
        return;
    }
    q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
    if (isfinite(q->z)) {
        if (integer(q)) {
            bound(q,h,zp,x);
        } else if (branch(q, *zp)) {
            push(h, q);
            return;
        }
    }
    free_node(q);
}

double intopt(int m, int n, double** a, double* b, double* c, double* x) {
    struct node_t* p = initial_node(m, n, a, b, c);
    struct Node* h = NULL;
    
    push(&h, p);
    
    double z = -INFINITY;
    p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
    if (integer(p) || !isfinite(p->z)) {
        z = p->z;
        if (integer(p)) {
            memcpy(x, p->x, (p->n) * sizeof(double));
        }

        free_node(p);
        if (h != NULL) {
            free_stack(&h);
        }
        
        return z;
    }
    branch(p,z);
    while (h != NULL) {
        p = pop(&h);

        succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
        succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);

        free_node(p);
        
    }
    free_stack(&h);
    if(z == -INFINITY) {
        return NAN;
    } else {
        return z;
    }
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

/*

int main(int argc, char** argv) {
    int m, n;
    double* c;
    double** a;
    double* b;
    int i, j;
    int* var; 

    scanf("%d %d", &m, &n);
    printf("m = %d, n = %d\n", m, n);
    double y = 0;
    double* x = calloc((n+m), sizeof(double)); 

   

    c = calloc(n, sizeof(double));
    
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

    print_matrix(m, n, c, a, b);


    double min_index = intopt(m, n, a, b, c, x);

    printf("Index of minimum b[i]: %lf\n", min_index);

    free(c);
    for (i = 0; i < m; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);

    return 0;
}*/
