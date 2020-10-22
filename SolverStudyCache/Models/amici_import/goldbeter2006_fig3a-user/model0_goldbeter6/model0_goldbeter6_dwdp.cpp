#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_goldbeter6(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*Q;
            break;
        case 1:
            dwdp[1] = -1.0*P*reaction_1_V/pow(1.0*P + reaction_1_Km, 2);
            break;
        case 2:
            dwdp[1] = 1.0*P/(1.0*P + reaction_1_Km);
            break;
        case 3:
            dwdp[2] = -reaction_2_V1*(-1.0*Q + 1)/pow(-1.0*Q + reaction_2_K1 + 1, 2);
            break;
        case 4:
            dwdp[2] = (-1.0*Q + 1)/(-1.0*Q + reaction_2_K1 + 1);
            break;
        case 5:
            dwdp[3] = -1.0*Q*R*reaction_3_V2/pow(1.0*Q + reaction_3_K2, 2);
            break;
        case 6:
            dwdp[3] = 1.0*Q*R/(1.0*Q + reaction_3_K2);
            break;
        case 7:
            dwdp[4] = -1.0*P*reaction_4_V3*(-1.0*R + 1)/pow(-1.0*R + reaction_4_k3 + 1, 2);
            break;
        case 8:
            dwdp[4] = 1.0*P*(-1.0*R + 1)/(-1.0*R + reaction_4_k3 + 1);
            break;
        case 9:
            dwdp[5] = 1.0*R/(1.0*R + reaction_5_Km);
            break;
        case 10:
            dwdp[5] = -1.0*R*reaction_5_V/pow(1.0*R + reaction_5_Km, 2);
            break;
    }
}