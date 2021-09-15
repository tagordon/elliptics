#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

//void bdj_arrays(double s[], double n[], double mc[], double *b, double *d, double *j, int *k);
//void fkp(double phi[], double n[], double mc[], double *f, double *e, double *p, int *k);
//void integrate_linear(double *s1, double *s2, double *r, double *b, double *f);
void f_burl(double *phi, double *m, double *e);
void e_burl(double *phi, double *m, double *e);
void p_burl(double *phi, double *n, double *m, double *e);

void kc_burl(double *k, double *e);
void ec_burl(double *k, double *e);
void pc_burl(double *k, double *p, double *e);
//void test_int(double phi1[], double phi2[], double b[], double r[], double *res, int *j);
void flux(double *c1, double *c2, double *rp, double *rm, double bp2[], double bm2[], double bpm2[], double *lc, int *j);
void Arc(double *c1, double *c2, double *phi1, double *phi2, double *r, double *b, double *res);
void F_lin(double *phi, double *r, double *b, double *res);
void F_const(double *phi, double *r, double *b, double *res);
void F_quad(double *phi, double *r, double *b, double *res);

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

double Flin(double phi, double r, double b, double res){
    F_lin(&phi, &r, &b, &res);
    return res;
}

double Fconst(double phi, double r, double b, double res){
    F_const(&phi, &r, &b, &res);
    return res;
}

double Fquad(double phi, double r, double b, double res){
    F_quad(&phi, &r, &b, &res);
    return res;
}

//void I(double *phi1, double *phi2, double *b, double *r, double *res, int j){
//    test_int(phi1, phi2, r, b, res, &j);
//}

double ellF(double phi, double m, double e){
    f_burl(&phi, &m, &e);
    return e;
}

double ellE(double phi, double m, double e){
    e_burl(&phi, &m, &e);
    return e;
}

double ellP(double phi, double n, double m, double e){
    p_burl(&phi, &n, &m, &e);
    return e;
}

double ellPC(double k, double p, double e){
    pc_burl(&k, &p, &e);
    return e;
}

double ellEC(double k, double e){
    ec_burl(&k, &e);
    return e;
}

double ellKC(double k, double e){
    kc_burl(&k, &e);
    return e;
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