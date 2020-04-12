#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_ma(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = k11;
    dwdx[1] = PKA*k2;
    dwdx[2] = k9;
    dwdx[3] = k1;
    dwdx[4] = k14;
    dwdx[5] = k5;
    dwdx[6] = PKA*k6;
    dwdx[7] = REGA*k8;
    dwdx[8] = ACA*k2;
    dwdx[9] = k4;
    dwdx[10] = ERK2*k6;
    dwdx[11] = incAMP*k10;
    dwdx[12] = ERK2*k8;
    dwdx[13] = k12;
    dwdx[14] = k13;
    dwdx[15] = REGA*k10;
    dwdx[16] = k3;
}