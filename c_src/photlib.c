#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void flux(double *c1, double *c2, double *rp, double *rm, double bp2[], double bm2[], double cth[], double sth[], double **lc, int *j);

void coords(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double xs[], double ys[], double zs[], double xp[], double yp[], double zp[], double xm[], double ym[], double zm[]);

void impacts(double t[], double *ms, double *t0p, double *ep, double *Pp, double *Op, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double x[], double y[], double xbc[], double ybc[], double bp2[], double bm2[], double bpm2[]);

void bp_bpm_theta(double t[], double *ms, double *t0p, double *ep, double *Pp, double *wp, double *ip, double *mp, \
    double *t0m, double *em, double *Pm, double *Om, double *wm, double *im, double *mm, int *j, \
    double bp[], double bpm[], double theta[]);

int main(){
    return 0;
}

void system_coords(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *xs, double *ys, double *zs, double *xp, double *yp, double *zp, double *xm, double *ym, double *zm){
    coords(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, xs, ys, zs, xp, yp, zp, xm, ym, zm);
}

void system_impacts(double *t, double ms, double t0p, double ep, double Pp, double Op, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *x, double *y, double *xbc, double*ybc, double *bp2, double *bm2, double *bpm2){
    impacts(t, &ms, &t0p, &ep, &Pp, &Op, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, x, y, xbc, ybc, bp2, bm2, bpm2);
}

void system_bp_bpm_theta(double *t, double ms, double t0p, double ep, double Pp, double wp, double ip, double mp, \
    double t0m, double em, double Pm, double Om, double wm, double im, double mm, int j, \
    double *bp, double *bpm, double *theta){
    bp_bpm_theta(t, &ms, &t0p, &ep, &Pp, &wp, &ip, &mp, &t0m, &em, &Pm, &Om, &wm, &im, &mm, &j, bp, bpm, theta);
}

void LC(double c1, double c2, double rp, double rm, double *bp, double *bpm, double *cth, double*sth, double **lc, int j){
    flux(&c1, &c2, &rp, &rm, bp, bpm, cth, sth, lc, &j);
}