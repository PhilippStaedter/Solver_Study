#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_aguda1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = aCYCEcdk2*pRBE2F;
            break;
        case 1:
            dwdp[5] = CYCEcdk2p27;
            break;
        case 2:
            dwdp[6] = CYCDcdk4*aCYCEcdk20;
            break;
        case 3:
            dwdp[7] = E2F;
            break;
        case 4:
            dwdp[8] = CYCDcdk4*p27;
            break;
        case 5:
            dwdp[0] = CYCDcdk4*pRBE2F;
            break;
        case 6:
            dwdp[0] = CYCDcdk4p27*pRBE2F;
            break;
        case 7:
            dwdp[20] = aCYCEcdk2*iCYCEcdk2;
            break;
        case 8:
            dwdp[9] = CYCDcdk4p27;
            break;
        case 9:
            dwdp[10] = pow(aCYCEcdk2, 2);
            break;
        case 10:
            dwdp[12] = p27;
            break;
        case 11:
            dwdp[13] = 1;
            break;
        case 12:
            dwdp[14] = aCYCEcdk20;
            break;
        case 13:
            dwdp[15] = 1.0/(k25a*pRB + 1);
            break;
        case 14:
            dwdp[15] = -k25*pRB/pow(k25a*pRB + 1, 2);
            break;
        case 15:
            dwdp[16] = 1.0/(aCYCEcdk20*k26a + 1);
            break;
        case 16:
            dwdp[16] = -aCYCEcdk20*k26/pow(aCYCEcdk20*k26a + 1, 2);
            break;
        case 17:
            dwdp[17] = 1;
            break;
        case 18:
            dwdp[18] = pRB;
            break;
        case 19:
            dwdp[19] = aCYCEcdk21;
            break;
        case 20:
            dwdp[22] = E2F;
            break;
        case 21:
            dwdp[22] = 1;
            break;
        case 22:
            dwdp[23] = 1;
            break;
        case 23:
            dwdp[25] = iCYCEcdk2;
            break;
        case 24:
            dwdp[26] = 1;
            break;
        case 25:
            dwdp[2] = 1;
            break;
        case 26:
            dwdp[3] = aCYCEcdk2*p27;
            break;
        case 27:
            dwdp[4] = aCYCEcdk2*p27;
            break;
        case 28:
            dwdp[11] = E2F*pRB;
            break;
        case 29:
            dwdp[21] = aCYCEcdk2;
            break;
        case 30:
            dwdp[24] = E2F;
            break;
        case 31:
            dwdp[1] = CYCDcdk4;
            break;
    }
}