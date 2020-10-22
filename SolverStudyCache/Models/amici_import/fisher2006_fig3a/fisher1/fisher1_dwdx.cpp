#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_fisher1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -2.6900000000000001e-13*NFAT_Cyt*k16;
    dwdx[1] = -2.6900000000000001e-13*k20;
    dwdx[2] = -2.6900000000000001e-13*k5;
    dwdx[3] = -2.6900000000000001e-13*NFAT_Pi_Cyt*k11;
    dwdx[4] = -1.13e-13*k20;
    dwdx[5] = 1.13e-13*NFAT_Nuc*k16;
    dwdx[6] = 1.13e-13*k6;
    dwdx[7] = -1.13e-13*NFAT_Pi_Nuc*k11;
    dwdx[8] = 8.070000000000001e-13*pow(Ca_Cyt, 2)*Inact_C_Cyt*k19;
    dwdx[9] = 2.6900000000000001e-13*k21;
    dwdx[10] = 3.3899999999999997e-13*pow(Ca_Nuc, 2)*Inact_C_Nuc*k19;
    dwdx[11] = -1.13e-13*k22;
    dwdx[12] = 2.6900000000000001e-13*pow(Ca_Cyt, 3)*k19;
    dwdx[13] = 2.6900000000000001e-13*k5;
    dwdx[14] = 1.13e-13*pow(Ca_Nuc, 3)*k19;
    dwdx[15] = -1.13e-13*k6;
    dwdx[16] = 2.6900000000000001e-13*k15;
    dwdx[17] = -2.6900000000000001e-13*k9;
    dwdx[18] = 2.6900000000000001e-13*k14;
    dwdx[19] = -1.13e-13*k15;
    dwdx[20] = 1.13e-13*k14;
    dwdx[21] = 1.13e-13*k10;
    dwdx[22] = -2.6900000000000001e-13*Act_C_Cyt*k16;
    dwdx[23] = -2.6900000000000001e-13*k2;
    dwdx[24] = -2.6900000000000001e-13*k17;
    dwdx[25] = -1.13e-13*k2;
    dwdx[26] = 1.13e-13*Act_C_Nuc*k16;
    dwdx[27] = 1.13e-13*k18;
    dwdx[28] = 2.6900000000000001e-13*k7;
    dwdx[29] = -2.6900000000000001e-13*k13;
    dwdx[30] = 2.6900000000000001e-13*k12;
    dwdx[31] = -1.13e-13*k8;
    dwdx[32] = -1.13e-13*k13;
    dwdx[33] = 1.13e-13*k12;
    dwdx[34] = 2.6900000000000001e-13*k3;
    dwdx[35] = 2.6900000000000001e-13*k1;
    dwdx[36] = -2.6900000000000001e-13*Act_C_Cyt*k11;
    dwdx[37] = 1.13e-13*k1;
    dwdx[38] = -1.13e-13*k4;
    dwdx[39] = -1.13e-13*Act_C_Nuc*k11;
}