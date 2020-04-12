#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_lee4(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*P*k1;
    dwdx[1] = 1.0*P21*k7;
    dwdx[2] = 1.0*P2*k7a;
    dwdx[3] = 1.0*M1*k3;
    dwdx[4] = 1.0*M*k3a;
    dwdx[5] = 1.0*P*k5;
    dwdx[6] = -1.0*j3a;
    dwdx[7] = 1.0*k4a;
    dwdx[8] = -1.0*j3;
    dwdx[9] = 1.0*k4;
    dwdx[10] = -1.0*j7a;
    dwdx[11] = 1.0*k8a;
    dwdx[12] = -1.0*j7;
    dwdx[13] = 1.0*k8;
    dwdx[14] = -1.0*j1;
    dwdx[15] = 1.0*kC1;
    dwdx[16] = 1.0*k2;
    dwdx[17] = 1.0*kC2;
    dwdx[18] = -1.0*j5;
    dwdx[19] = 1.0*k6;
    dwdx[20] = 1.0*k9;
    dwdx[21] = 1.0*E*k3a;
    dwdx[22] = 1.0*E*k3;
    dwdx[23] = 1.0*E*k1;
    dwdx[24] = 1.0*E*k5;
    dwdx[25] = 1.0*k10;
    dwdx[26] = 1.0*E*k7a;
    dwdx[27] = 1.0*E*k7;
}