#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model2_balagadde1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*C1*kc1*(-C1 - C2)/pow(Cm, 2);
            dwdp[2] = -1.0*C2*kc2*(-C1 - C2)/pow(Cm, 2);
            break;
        case 1:
            dwdp[1] = 1.0*C1;
            dwdp[3] = 1.0*C2;
            dwdp[5] = 1.0*A1;
            dwdp[7] = 1.0*A2;
            break;
        case 2:
            dwdp[1] = 1.0*C1*(-K1*(pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.5)/pow(pow(A2, 2) + K1, 2) + (pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.5)/(pow(A2, 2) + K1));
            break;
        case 3:
            dwdp[3] = -1.0*pow(A1, 2)*C2*d2/pow(pow(A1, 2) + K2, 2);
            break;
        case 4:
            dwdp[3] = 1.0*pow(A1, 2)*C2/(pow(A1, 2) + K2);
            break;
        case 5:
            dwdp[5] = 1.0*A1;
            break;
        case 6:
            dwdp[7] = 1.0*A2;
            break;
        case 7:
            dwdp[4] = 1.0*C1;
            break;
        case 8:
            dwdp[0] = 1.0*C1*(1 - (C1 + C2)/Cm);
            break;
        case 9:
            dwdp[2] = 1.0*C2*(1 - (C1 + C2)/Cm);
            break;
    }
}