#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846

void bjd(double s0, double n, double mc, double b, double d, double j);
void serj(double y, double n, double m);
void uatan(double t, double h);

void FJ(double * res, double phi, double n, double m){
    
    double x, y;
    double d, u, J;
    double h, T;
    double s, c;
    int N;
    
    if (phi > PI/2.){
        res[0] = NAN;
        res[1] = NAN;
        return;
    }
    
    s = sin(phi);
    c = cos(phi);
    y = pow(s, 2.);
    x = pow(c, 2.);
    h = n * (1 - n) * (n - m);
    
    N = 0;
    while (y > 0.5){
        d = sqrt(1 - m + m*x);
        x = (sqrt(x) + d) / (1 + d);
        y = 1 - x;
        ++N;
    }
    while (y > 0.001){
        y = y / ((1 + sqrt(1 - y)) * (1 + sqrt(1 - m*y)));
        ++N;
    }
    
    u = 1 + y * (1 + m)/6. 
        + pow(y, 2.) * (3 + 2*m + 3*pow(m, 2.))/40.
        + pow(y, 3.) * (5 + 3*m + 3*pow(m, 2.) + 5*pow(m, 3.))/112.
        + pow(y, 4.) * (35 + 20*m + 18*pow(m, 2.) + 20*pow(m, 3.) + 35*pow(m, 4.))/1152. 
        + pow(y, 5.) * (63 + 35*m + 30*pow(m, 2.) + 30*pow(m, 3.) + 35*pow(m, 4.) + 63*pow(m, 5.))/2816. 
        + pow(y, 6.) * (231 + 126*m + 105*pow(m, 2.) + 100*pow(m, 3.) 
                        + 105*pow(m, 4.) + 126*pow(m, 5.) + 231*pow(m, 6.))/13312.;
    
    u = u * sqrt(y);
    res[0] = (1 << N) * u * ((phi > 0) - (phi < 0));
    
    J = pow(u, 3.) * (1 / 3.) - pow(u, 5.) * (1 + m - 3 * n)/15.
        + pow(u, 7.) *  (2 + 13*m + 2*pow(m, 2.) - n*(30+30*m) + 45*pow(n, 2)) /315.
        - pow(u, 9.) *  (1 + 30*m + 30*pow(m, 2.) + pow(m, 3.) 
                         - (63 + 252*m + 63*pow(m, 2))*n + (315 + 315*m)*pow(n, 2.) 
                         - 315*pow(n, 3.))/2835.
        + pow(u, 11.) *  (2 + 251*m + 876*pow(m, 2.) + 251*pow(m, 3.) 
                          + 2*pow(m, 4.) - (510 + 5850*m + 5850*pow(m, 2.) + 510*pow(m, 3.))*n 
                          + (6615 + 21735*m + 6615*pow(m, 2.))*pow(n, 2.) 
                          - (18900 + 18900*m)*pow(n, 3.) 
                          + 14175*pow(n, 4.))/155925.
        - pow(u, 13.) * (2 + 1018*m + 9902*pow(m, 2.) + 9902*pow(m, 3.) + 1018*pow(m, 4.) + 2*pow(m, 5.)
                        - (2046 + 59268*m + 158103*pow(m, 2.) + 59268*pow(m, 3.) + 2046*pow(m, 4.))*n
                        + (63360 + 497475*m + 497475*pow(m, 2.) + 63360*pow(m, 3.) )*pow(n, 2.)
                        - (395010 + 1164240*m + 395010*pow(m, 2.))*pow(n, 3.)
                        + (779625 + 779625*m)*pow(n, 4.) - 467775*pow(n, 5.))/6081075.;
    
    if (h > 0){
        T = atan2(sqrt(h)*s, c)/sqrt(h);
    }
    else if (h == 0){
        T = s / c;
    }
    else{
        h = -h;
        T = atanh(sqrt(h)*s/c)/sqrt(h);
    }
    
    //for (i=0;i<N;i++){
    //    J = 2 * J + T;
    //}
    res[1] = J;
}
