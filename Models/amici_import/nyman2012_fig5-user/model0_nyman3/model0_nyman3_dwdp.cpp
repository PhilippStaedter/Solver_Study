#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_nyman3(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = IR*k1a;
            break;
        case 1:
            dwdp[0] = IR*ins;
            break;
        case 2:
            dwdp[0] = IR;
            break;
        case 3:
            dwdp[1] = IRins;
            break;
        case 4:
            dwdp[2] = IRins;
            break;
        case 5:
            dwdp[3] = IRp;
            break;
        case 6:
            dwdp[4] = IRiP;
            break;
        case 7:
            dwdp[4] = IRiP*Xp/(Xp + 1);
            break;
        case 8:
            dwdp[5] = IRp;
            break;
        case 9:
            dwdp[6] = IRi;
            break;
        case 10:
            dwdp[7] = IRS*(IRiP*k22 + IRp)/(Xp*km23 + 1);
            break;
        case 11:
            dwdp[7] = IRS*IRiP*k21/(Xp*km23 + 1);
            break;
        case 12:
            dwdp[8] = IRSiP*X;
            break;
        case 13:
            dwdp[9] = IRSiP;
            break;
        case 14:
            dwdp[7] = -IRS*Xp*k21*(IRiP*k22 + IRp)/pow(Xp*km23 + 1, 2);
            break;
        case 15:
            dwdp[10] = Xp;
            break;
    }
}