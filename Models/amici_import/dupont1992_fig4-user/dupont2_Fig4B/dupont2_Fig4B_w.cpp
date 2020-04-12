#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_dupont2_Fig4B(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*vp*(-W_star/(K2 + W_star) + pow(Z, q)*vMK*(-W_star + 1)/(vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)))/Wt;
    w[1] = 1.0*v0;
    w[2] = 1.0*v1_beta;
    w[3] = 1.0*Vm2*pow(Z, n)/(pow(Kp, n) + pow(Z, n));
    w[4] = 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)/((pow(K_A, amici_p) + pow(Z, amici_p))*(pow(Kr, m) + pow(Y, m)));
    w[5] = 1.0*Z*amici_k;
    w[6] = 1.0*Y*kf;
}