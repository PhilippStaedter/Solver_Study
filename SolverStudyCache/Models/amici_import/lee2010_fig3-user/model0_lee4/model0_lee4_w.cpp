#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_lee4(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*E*P*k1 - 1.0*E_P_1*j1;
    w[1] = 1.0*P2*k10;
    w[2] = 1.0*E*P21*k7 - 1.0*E_P21*j7;
    w[3] = 1.0*E*P2*k7a - 1.0*E_P2*j7a;
    w[4] = 1.0*E_P21*k8;
    w[5] = 1.0*E_P2*k8a;
    w[6] = 1.0*E_P_1*kC1;
    w[7] = 1.0*E_P_2*kC2;
    w[8] = 1.0*E_P_1*k2;
    w[9] = 1.0*M*k9;
    w[10] = 1.0*E*M1*k3 - 1.0*E_M1*j3;
    w[11] = 1.0*E*M*k3a - 1.0*E_M*j3a;
    w[12] = 1.0*E_M1*k4;
    w[13] = 1.0*E_M*k4a;
    w[14] = 1.0*E*P*k5 - 1.0*E_P_2*j5;
    w[15] = 1.0*E_P_2*k6;
}