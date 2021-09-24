#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void flux(double *c1, double *c2, double *rp, double *rm, double bp2[], double bm2[], double bpm2[], double *lc, int *j);

void coords(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double xs[], double ys[], double zs[], double xp[], double yp[], double zp[], double xm[], double ym[], double zm[]);

int main(){
    return 0;
}

void system_coords(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *xs, double *ys, double *zs, double *xp, double *yp, double *zp, double *xm, double *ym, double *zm){
    coords(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, xs, ys, zs, xp, yp, zp, xm, ym, zm);
}

void LC(double c1, double c2, double rp, double rm, double *bp2, double *bm2, double *bpm2, double *lc, int j){
    flux(&c1, &c2, &rp, &rm, bp2, bm2, bpm2, lc, &j);
}