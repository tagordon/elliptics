#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

//void bdj_arrays(double s[], double n[], double mc[], double *b, double *d, double *j, int *k);
//void fkp(double phi[], double n[], double mc[], double *f, double *e, double *p, int *k);
//void integrate_linear(double *s1, double *s2, double *r, double *b, double *f);
//void f_burl(double phi[], double m[], double *e, int *j);
//void e_burl(double phi[], double m[], double *e, int *j);
//void p_burl(double phi[], double n[], double m[], double *e, int *j);
//void test_int(double phi1[], double phi2[], double b[], double r[], double *res, int *j);
void flux(double *c1, double *c2, double *rp, double *rm, double bp2[], double bm2[], double bpm2[], double *lc, int *j);
void Arc(double *c1, double *c2, double *phi1, double *phi2, double *r, double *b, double *res);

int main(){
    return 0;
}

void LC(double c1, double c2, double rp, double rm, double *bp2, double *bm2, double *bpm2, double *lc, int j){
    flux(&c1, &c2, &rp, &rm, bp2, bm2, bpm2, lc, &j);
}

double Area(double c1, double c2, double r, double b, double e, double res){
    double phi1 = -PI * e;
    double phi2 = PI * e;
    Arc(&c1, &c2, &phi1, &phi2, &r, &b, &res);
    return res;
}

//void I(double *phi1, double *phi2, double *b, double *r, double *res, int j){
//    test_int(phi1, phi2, r, b, res, &j);
//}

//void F(double *phi, double *m, double *e, int j){
//    f_burl(phi, m, e, &j);
//}

//void E(double *phi, double *m, double *e, int j){
//    e_burl(phi, m, e, &j);
//}

//void P(double *phi, double *n, double *m, double *e, int j){
//    p_burl(phi, n, m, e, &j);
//}

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