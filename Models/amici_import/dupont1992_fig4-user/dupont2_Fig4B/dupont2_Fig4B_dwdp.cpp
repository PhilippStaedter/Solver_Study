#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_dupont2_Fig4B(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*pow(Z, q)*vMK*(-W_star + 1)/(Wt*(pow(Ka, q) + pow(Z, q))*pow(K1 - W_star + 1, 2));
            break;
        case 1:
            dwdp[0] = 1.0*W_star*vp/(Wt*pow(K2 + W_star, 2));
            break;
        case 2:
            dwdp[4] = -1.0*pow(K_A, amici_p)*Vm3*pow(Y, m)*pow(Z, amici_p)*amici_p/(K_A*pow(pow(K_A, amici_p) + pow(Z, amici_p), 2)*(pow(Kr, m) + pow(Y, m)));
            break;
        case 3:
            dwdp[0] = -1.0*pow(Ka, q)*pow(Z, q)*q*vMK*(-W_star + 1)/(Ka*Wt*pow(pow(Ka, q) + pow(Z, q), 2)*(K1 - W_star + 1));
            break;
        case 4:
            dwdp[3] = -1.0*pow(Kp, n)*Vm2*pow(Z, n)*n/(Kp*pow(pow(Kp, n) + pow(Z, n), 2));
            break;
        case 5:
            dwdp[4] = -1.0*pow(Kr, m)*Vm3*pow(Y, m)*pow(Z, amici_p)*m/(Kr*(pow(K_A, amici_p) + pow(Z, amici_p))*pow(pow(Kr, m) + pow(Y, m), 2));
            break;
        case 6:
            dwdp[3] = 1.0*pow(Z, n)/(pow(Kp, n) + pow(Z, n));
            break;
        case 7:
            dwdp[4] = 1.0*pow(Y, m)*pow(Z, amici_p)/((pow(K_A, amici_p) + pow(Z, amici_p))*(pow(Kr, m) + pow(Y, m)));
            break;
        case 8:
            dwdp[5] = 1.0*Z;
            break;
        case 9:
            dwdp[6] = 1.0*Y;
            break;
        case 10:
            dwdp[4] = 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)*log(Y)/((pow(K_A, amici_p) + pow(Z, amici_p))*(pow(Kr, m) + pow(Y, m))) + 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)*(-pow(Kr, m)*log(Kr) - pow(Y, m)*log(Y))/((pow(K_A, amici_p) + pow(Z, amici_p))*pow(pow(Kr, m) + pow(Y, m), 2));
            break;
        case 11:
            dwdp[3] = 1.0*Vm2*pow(Z, n)*log(Z)/(pow(Kp, n) + pow(Z, n)) + 1.0*Vm2*pow(Z, n)*(-pow(Kp, n)*log(Kp) - pow(Z, n)*log(Z))/pow(pow(Kp, n) + pow(Z, n), 2);
            break;
        case 12:
            dwdp[4] = 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)*log(Z)/((pow(K_A, amici_p) + pow(Z, amici_p))*(pow(Kr, m) + pow(Y, m))) + 1.0*Vm3*pow(Y, m)*pow(Z, amici_p)*(-pow(K_A, amici_p)*log(K_A) - pow(Z, amici_p)*log(Z))/(pow(pow(K_A, amici_p) + pow(Z, amici_p), 2)*(pow(Kr, m) + pow(Y, m)));
            break;
        case 13:
            dwdp[0] = 1.0*vp*(pow(Z, q)*vMK*(-W_star + 1)*log(Z)/(vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)) + pow(Z, q)*vMK*(-W_star + 1)*(-pow(Ka, q)*log(Ka) - pow(Z, q)*log(Z))/(vp*pow(pow(Ka, q) + pow(Z, q), 2)*(K1 - W_star + 1)))/Wt;
            break;
        case 14:
            dwdp[1] = 1.0;
            break;
        case 15:
            dwdp[2] = 1.0;
            break;
        case 16:
            dwdp[0] = 1.0*pow(Z, q)*(-W_star + 1)/(Wt*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1));
            break;
        case 17:
            dwdp[0] = -1.0*pow(Z, q)*vMK*(-W_star + 1)/(Wt*vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)) + 1.0*(-W_star/(K2 + W_star) + pow(Z, q)*vMK*(-W_star + 1)/(vp*(pow(Ka, q) + pow(Z, q))*(K1 - W_star + 1)))/Wt;
            break;
    }
}