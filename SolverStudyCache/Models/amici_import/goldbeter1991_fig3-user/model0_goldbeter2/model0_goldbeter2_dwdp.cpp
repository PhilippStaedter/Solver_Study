#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_goldbeter2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[3] = -Cvar*VM1*cell*(-Mvar + 1)/((Cvar + Kc)*pow(K1 - Mvar + 1, 2));
            break;
        case 1:
            dwdp[4] = -Mvar*V2*cell/pow(K2 + Mvar, 2);
            break;
        case 2:
            dwdp[5] = -Mvar*VM3*cell*(-X + 1)/pow(K3 - X + 1, 2);
            break;
        case 3:
            dwdp[6] = -V4*X*cell/pow(K4 + X, 2);
            break;
        case 4:
            dwdp[3] = -Cvar*VM1*cell*(-Mvar + 1)/(pow(Cvar + Kc, 2)*(K1 - Mvar + 1));
            break;
        case 5:
            dwdp[2] = -Cvar*X*cell*vd/pow(Cvar + Kd, 2);
            break;
        case 6:
            dwdp[4] = Mvar*cell/(K2 + Mvar);
            break;
        case 7:
            dwdp[6] = X*cell/(K4 + X);
            break;
        case 8:
            dwdp[3] = Cvar*cell*(-Mvar + 1)/((Cvar + Kc)*(K1 - Mvar + 1));
            break;
        case 9:
            dwdp[5] = Mvar*cell*(-X + 1)/(K3 - X + 1);
            break;
        case 10:
            dwdp[0] = vi;
            dwdp[1] = Cvar*kd;
            dwdp[2] = Cvar*X*vd/(Cvar + Kd);
            dwdp[3] = Cvar*VM1*(-Mvar + 1)/((Cvar + Kc)*(K1 - Mvar + 1));
            dwdp[4] = Mvar*V2/(K2 + Mvar);
            dwdp[5] = Mvar*VM3*(-X + 1)/(K3 - X + 1);
            dwdp[6] = V4*X/(K4 + X);
            break;
        case 11:
            dwdp[1] = Cvar*cell;
            break;
        case 12:
            dwdp[2] = Cvar*X*cell/(Cvar + Kd);
            break;
        case 13:
            dwdp[0] = cell;
            break;
    }
}