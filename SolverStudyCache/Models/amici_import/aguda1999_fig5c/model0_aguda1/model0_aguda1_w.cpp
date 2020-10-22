#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_aguda1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = CYCDcdk4*k1a*pRBE2F + CYCDcdk4p27*k1aa*pRBE2F + aCYCEcdk2*k1*pRBE2F;
    w[1] = CYCDcdk4*kmin6;
    w[2] = k7;
    w[3] = aCYCEcdk2*k8*p27;
    w[4] = aCYCEcdk2*k9*p27;
    w[5] = CYCEcdk2p27*k10;
    w[6] = CYCDcdk4*aCYCEcdk20*k17;
    w[7] = E2F*k18;
    w[8] = CYCDcdk4*k19*p27;
    w[9] = CYCDcdk4p27*k20;
    w[10] = pow(aCYCEcdk2, 2)*k21;
    w[11] = E2F*kmin1*pRB;
    w[12] = k22*p27;
    w[13] = k23;
    w[14] = aCYCEcdk20*k24;
    w[15] = k25/(k25a*pRB + 1);
    w[16] = k26/(aCYCEcdk20*k26a + 1);
    w[17] = k27;
    w[18] = k28*pRB;
    w[19] = aCYCEcdk21*k29;
    w[20] = aCYCEcdk2*iCYCEcdk2*k2;
    w[21] = aCYCEcdk2*kmin2;
    w[22] = E2F*k3 + k3a;
    w[23] = k4;
    w[24] = E2F*kmin4;
    w[25] = iCYCEcdk2*k5;
    w[26] = k6;
}