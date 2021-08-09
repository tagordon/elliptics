#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

//void bdj_arrays(double s[], double n[], double mc[], double *b, double *d, double *j, int *k);
//void fkp(double phi[], double n[], double mc[], double *f, double *e, double *p, int *k);
//void integrate_linear(double *s1, double *s2, double *r, double *b, double *f);
void f_burl(double phi[], double m[], double *e, int *j);
void e_burl(double phi[], double m[], double *e, int *j);
void p_burl(double phi[], double n[], double m[], double *e, int *j);

int main(){
    return 0;
}

void F(double *phi, double *m, double *e, int j){
    f_burl(phi, m, e, &j);
}

void E(double *phi, double *m, double *e, int j){
    e_burl(phi, m, e, &j);
}

void P(double *phi, double *n, double *m, double *e, int j){
    p_burl(phi, n, m, e, &j);
}

//double integrate_along_curve(double s1, double s2, double r, double b){
//    double f;
//    integrate_linear(&s1, &s2, &r, &b, &f);
//    return f;
//}

//void BDJ_arrays(double *b, double *d, double *j, double * s, double *n, double *m, int k){
//    bdj_arrays(s, n, m, b, d, j, &k);
//}

//void FKP_arrays(double *f, double *e, double *p, double *phi, double *n, double *m, int k){
//    fkp(phi, n, m, f, e, p, &k);
//}