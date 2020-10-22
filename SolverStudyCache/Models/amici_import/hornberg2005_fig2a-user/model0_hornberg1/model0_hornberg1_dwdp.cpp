#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_hornberg1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[7] = Inh*Vm8*x3p/(pow(Ki8, 2)*Km8*pow(Inh/Ki8 + 1 + x3p/Km8, 2));
            break;
        case 1:
            dwdp[0] = -R*Vm1/pow(Km1 + R, 2);
            break;
        case 2:
            dwdp[1] = -Rin*Vm2/pow(Km2 + Rin, 2);
            break;
        case 3:
            dwdp[2] = -R*k3*x1/pow(Km3 + x1, 2);
            break;
        case 4:
            dwdp[3] = -Vm4*x1p/pow(Km4 + x1p, 2);
            break;
        case 5:
            dwdp[4] = -k5*x1p*x2/pow(Km5 + x2, 2);
            break;
        case 6:
            dwdp[5] = -Vm6*x2p/pow(Km6 + x2p, 2);
            break;
        case 7:
            dwdp[6] = -k7*x2p*x3/pow(Km7 + x3, 2);
            break;
        case 8:
            dwdp[7] = -Vm8*x3p/(pow(Km8, 2)*(Inh/Ki8 + 1 + x3p/Km8)) + Vm8*pow(x3p, 2)/(pow(Km8, 3)*pow(Inh/Ki8 + 1 + x3p/Km8, 2));
            break;
        case 9:
            dwdp[0] = R/(Km1 + R);
            break;
        case 10:
            dwdp[1] = Rin/(Km2 + Rin);
            break;
        case 11:
            dwdp[3] = x1p/(Km4 + x1p);
            break;
        case 12:
            dwdp[5] = x2p/(Km6 + x2p);
            break;
        case 13:
            dwdp[7] = x3p/(Km8*(Inh/Ki8 + 1 + x3p/Km8));
            break;
        case 14:
            dwdp[2] = R*x1/(Km3 + x1);
            break;
        case 15:
            dwdp[4] = x1p*x2/(Km5 + x2);
            break;
        case 16:
            dwdp[6] = x2p*x3/(Km7 + x3);
            break;
    }
}