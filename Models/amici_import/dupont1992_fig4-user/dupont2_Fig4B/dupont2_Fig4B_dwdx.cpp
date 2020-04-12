#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_dupont2_Fig4B(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*vp*(W_star/pow(K2 + W_star, 2) + pow(Z, q)*vMK*(-W_star + 1)/(vp*(pow(Ka, q) + pow(Z, q))*pow(K1 - W_star + 1, 2)) - pow(Z, q)*vMK/(vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)) - 1/(K2 + W_star))/Wt;
    dwdx[1] = -1.0*vp*(-W_star/(K2 + W_star) + pow(Z, q)*vMK*(-W_star + 1)/(vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)))/pow(Wt, 2);
    dwdx[2] = -1.0*Vm3*pow(Y, 2*m)*pow(Z, amici_p)*m/(Y*(pow(K_A, amici_p) + pow(Z, amici_p))*pow(pow(Kr, m) + pow(Y, m), 2)) + 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)*m/(Y*(pow(K_A, amici_p) + pow(Z, amici_p))*(pow(Kr, m) + pow(Y, m)));
    dwdx[3] = 1.0*kf;
    dwdx[4] = 1.0*vp*(-pow(Z, 2*q)*q*vMK*(-W_star + 1)/(Z*vp*pow(pow(Ka, q) + pow(Z, q), 2)*(K1 - W_star + 1)) + pow(Z, q)*q*vMK*(-W_star + 1)/(Z*vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)))/Wt;
    dwdx[5] = -1.0*Vm2*pow(Z, 2*n)*n/(Z*pow(pow(Kp, n) + pow(Z, n), 2)) + 1.0*Vm2*pow(Z, n)*n/(Z*(pow(Kp, n) + pow(Z, n)));
    dwdx[6] = -1.0*Vm3*pow(Y, m)*pow(Z, 2*amici_p)*amici_p/(Z*pow(pow(K_A, amici_p) + pow(Z, amici_p), 2)*(pow(Kr, m) + pow(Y, m))) + 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)*amici_p/(Z*(pow(K_A, amici_p) + pow(Z, amici_p))*(pow(Kr, m) + pow(Y, m)));
    dwdx[7] = 1.0*amici_k;
}