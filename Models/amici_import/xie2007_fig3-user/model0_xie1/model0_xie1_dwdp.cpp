#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_xie1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[17] = 1.0*CLK*CYC;
            break;
        case 1:
            dwdp[33] = 1.0*CC*(-prcpdp + 1);
            break;
        case 2:
            dwdp[28] = 1.0*CC*(-prcper + 1);
            break;
        case 3:
            dwdp[0] = 1.0*CC*PT;
            break;
        case 4:
            dwdp[38] = 1.0*CC*(-prct + 1);
            break;
        case 5:
            dwdp[31] = 1.0*CC*(-prcv + 1);
            break;
        case 6:
            dwdp[36] = 1.0*PDP*(-prpc - prvc + 1);
            break;
        case 7:
            dwdp[20] = 1.0*PER*TIM;
            break;
        case 8:
            dwdp[34] = 1.0*VRI*(-prpc - prvc + 1);
            break;
        case 9:
            dwdp[4] = 1.0*CC;
            break;
        case 10:
            dwdp[7] = 1.0*CCPT;
            break;
        case 11:
            dwdp[18] = 1.0*CLK;
            break;
        case 12:
            dwdp[40] = 1.0*clkm;
            break;
        case 13:
            dwdp[14] = 1.0*PDP;
            break;
        case 14:
            dwdp[12] = 1.0*pdpm;
            break;
        case 15:
            dwdp[22] = 1.0*PER;
            break;
        case 16:
            dwdp[2] = 1.0*perm;
            break;
        case 17:
            dwdp[15] = 1.0*PT;
            break;
        case 18:
            dwdp[26] = 1.0*TIM;
            break;
        case 19:
            dwdp[24] = 1.0*timm;
            break;
        case 20:
            dwdp[10] = 1.0*VRI;
            break;
        case 21:
            dwdp[8] = 1.0*vrim;
            break;
        case 22:
            dwdp[11] = 1.0*pdpp*(-tcccpdpp*pow(-prcpdp + 1, npdp)*log(-prcpdp + 1) + tcdvpmt*pow(-prcpdp + 1, npdp)*log(-prcpdp + 1));
            break;
        case 23:
            dwdp[5] = 1.0*perp*(-tcccperp*pow(-prcper + 1, npt)*log(-prcper + 1) + tcdvpmt*pow(-prcper + 1, npt)*log(-prcper + 1));
            dwdp[23] = 1.0*timp*(-tccctimp*pow(-prct + 1, npt)*log(-prct + 1) + tcdvpmt*pow(-prct + 1, npt)*log(-prct + 1));
            break;
        case 24:
            dwdp[6] = 1.0*vrip*(-tcccvrip*pow(-prcv + 1, nvri)*log(-prcv + 1) + tcdvpmt*pow(-prcv + 1, nvri)*log(-prcv + 1));
            break;
        case 25:
            dwdp[11] = 1.0*pdpp*(-pow(-prcpdp + 1, npdp) + 1);
            break;
        case 26:
            dwdp[5] = 1.0*perp*(-pow(-prcper + 1, npt) + 1);
            break;
        case 27:
            dwdp[23] = 1.0*timp*(-pow(-prct + 1, npt) + 1);
            break;
        case 28:
            dwdp[6] = 1.0*vrip*(-pow(-prcv + 1, nvri) + 1);
            break;
        case 29:
            dwdp[16] = 1.0*clkp*(-prpc - prvc + 1);
            break;
        case 30:
            dwdp[5] = 1.0*perp*pow(-prcper + 1, npt);
            dwdp[6] = 1.0*vrip*pow(-prcv + 1, nvri);
            dwdp[11] = 1.0*pdpp*pow(-prcpdp + 1, npdp);
            dwdp[23] = 1.0*timp*pow(-prct + 1, npt);
            break;
        case 31:
            dwdp[16] = 1.0*clkp*prpc;
            break;
        case 32:
            dwdp[16] = 1.0*clkp*prvc;
            break;
        case 33:
            dwdp[1] = 1.0*clkm;
            break;
        case 34:
            dwdp[13] = 1.0*pdpm;
            break;
        case 35:
            dwdp[3] = 1.0*perm;
            break;
        case 36:
            dwdp[25] = 1.0*timm;
            break;
        case 37:
            dwdp[9] = 1.0*vrim;
            break;
        case 38:
            dwdp[19] = 1.0*CC;
            break;
        case 39:
            dwdp[32] = 1.0*prcpdp;
            break;
        case 40:
            dwdp[29] = 1.0*prcper;
            break;
        case 41:
            dwdp[27] = 1.0*CCPT;
            break;
        case 42:
            dwdp[39] = 1.0*prct;
            break;
        case 43:
            dwdp[30] = 1.0*prcv;
            break;
        case 44:
            dwdp[37] = 1.0*prpc;
            break;
        case 45:
            dwdp[21] = 1.0*PT;
            break;
        case 46:
            dwdp[35] = 1.0*prvc;
            break;
    }
}