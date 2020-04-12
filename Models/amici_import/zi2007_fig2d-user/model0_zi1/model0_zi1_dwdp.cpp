#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_zi1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[13] = 0.0010499999999999999*LRC_EE;
            break;
        case 1:
            dwdp[0] = 0.0010499999999999999*T1R_EE;
            break;
        case 2:
            dwdp[6] = 0.0010499999999999999*T2R_EE;
            break;
        case 3:
            dwdp[16] = 0.00035*Smads_Complex_n;
            break;
        case 4:
            dwdp[18] = 0.00035*Smad2n;
            break;
        case 5:
            dwdp[20] = 0.00035*Smad4n;
            break;
        case 6:
            dwdp[10] = 0.0010499999999999999*Smad2c;
            break;
        case 7:
            dwdp[19] = 0.0010499999999999999*Smad4c;
            break;
        case 8:
            dwdp[15] = 0.0010499999999999999*Smads_Complex_c;
            break;
        case 9:
            dwdp[17] = 0.0010499999999999999*LRC_Cave*Smads_Complex_n;
            break;
        case 10:
            dwdp[7] = 0.0010499999999999999*T1R_Surf*T2R_Surf*TGF_beta;
            break;
        case 11:
            dwdp[14] = 0.0010499999999999999*LRC_EE*Smad2c*Smad4c;
            break;
        case 12:
            dwdp[2] = 0.0010499999999999999*T2R_Surf;
            dwdp[8] = 0.0010499999999999999*LRC_Surf;
            dwdp[22] = 0.0010499999999999999*T1R_Surf;
            break;
        case 13:
            dwdp[4] = 0.0010499999999999999*T2R_Surf;
            dwdp[11] = 0.0010499999999999999*LRC_Surf;
            dwdp[24] = 0.0010499999999999999*T1R_Surf;
            break;
        case 14:
            dwdp[3] = 0.0010499999999999999*T2R_Cave;
            dwdp[9] = 0.0010499999999999999*LRC_Cave;
            dwdp[23] = 0.0010499999999999999*T1R_Cave;
            break;
        case 15:
            dwdp[5] = 0.0010499999999999999*T2R_EE;
            dwdp[12] = 0.0010499999999999999*LRC_EE;
            dwdp[25] = 0.0010499999999999999*T1R_EE;
            break;
        case 16:
            dwdp[21] = 0.0010499999999999999;
            break;
        case 17:
            dwdp[1] = 0.0010499999999999999;
            break;
    }
}