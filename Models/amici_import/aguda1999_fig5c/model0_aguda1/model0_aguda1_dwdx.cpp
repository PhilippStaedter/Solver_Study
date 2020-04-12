#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_aguda1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = k1a*pRBE2F;
    dwdx[1] = kmin6;
    dwdx[2] = aCYCEcdk20*k17;
    dwdx[3] = k19*p27;
    dwdx[4] = k1aa*pRBE2F;
    dwdx[5] = k20;
    dwdx[6] = k10;
    dwdx[7] = k18;
    dwdx[8] = kmin1*pRB;
    dwdx[9] = k3;
    dwdx[10] = kmin4;
    dwdx[11] = k1*pRBE2F;
    dwdx[12] = k8*p27;
    dwdx[13] = k9*p27;
    dwdx[14] = 2*aCYCEcdk2*k21;
    dwdx[15] = iCYCEcdk2*k2;
    dwdx[16] = kmin2;
    dwdx[17] = CYCDcdk4*k17;
    dwdx[18] = k24;
    dwdx[19] = -k26*k26a/pow(aCYCEcdk20*k26a + 1, 2);
    dwdx[20] = k29;
    dwdx[21] = aCYCEcdk2*k2;
    dwdx[22] = k5;
    dwdx[23] = aCYCEcdk2*k8;
    dwdx[24] = aCYCEcdk2*k9;
    dwdx[25] = CYCDcdk4*k19;
    dwdx[26] = k22;
    dwdx[27] = E2F*kmin1;
    dwdx[28] = -k25*k25a/pow(k25a*pRB + 1, 2);
    dwdx[29] = k28;
    dwdx[30] = CYCDcdk4*k1a + CYCDcdk4p27*k1aa + aCYCEcdk2*k1;
}