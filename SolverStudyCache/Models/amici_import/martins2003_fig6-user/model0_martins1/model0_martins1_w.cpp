#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_martins1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = DFG*v1_k1;
    w[1] = E1*v10_k10;
    w[2] = E1*v11_k11;
    w[3] = Man*v12_k12;
    w[4] = Glu*v13_k13;
    w[5] = Cn*Gly*v14_k14;
    w[6] = Cn*v15_k15;
    w[7] = E2*v16_k16;
    w[8] = DFG*v2_k2;
    w[9] = DFG*v3_k3;
    w[10] = E1*v4_k4;
    w[11] = _3DG*v5_k5;
    w[12] = _3DG*v6_k6;
    w[13] = E2*v7_k7;
    w[14] = _1DG*v8_k8;
    w[15] = _1DG*v9_k9;
}