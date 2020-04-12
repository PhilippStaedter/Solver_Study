#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_martins1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = Gly*v14_k14;
    dwdx[1] = v15_k15;
    dwdx[2] = v1_k1;
    dwdx[3] = v2_k2;
    dwdx[4] = v3_k3;
    dwdx[5] = v10_k10;
    dwdx[6] = v11_k11;
    dwdx[7] = v4_k4;
    dwdx[8] = v16_k16;
    dwdx[9] = v7_k7;
    dwdx[10] = v13_k13;
    dwdx[11] = Cn*v14_k14;
    dwdx[12] = v12_k12;
    dwdx[13] = v8_k8;
    dwdx[14] = v9_k9;
    dwdx[15] = v5_k5;
    dwdx[16] = v6_k6;
}