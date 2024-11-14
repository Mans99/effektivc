#include<stdio.h>
#include<stdlib.h>

struct simplex_t {
    int m;
    int n;
    int* var;
    double ** A;
    double* x;
    double* c;
    double y;

};



