#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_lee4(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*E_P_1;
            break;
        case 1:
            dwdp[10] = -1.0*E_M1;
            break;
        case 2:
            dwdp[11] = -1.0*E_M;
            break;
        case 3:
            dwdp[14] = -1.0*E_P_2;
            break;
        case 4:
            dwdp[2] = -1.0*E_P21;
            break;
        case 5:
            dwdp[3] = -1.0*E_P2;
            break;
        case 6:
            dwdp[0] = 1.0*E*P;
            break;
        case 7:
            dwdp[1] = 1.0*P2;
            break;
        case 8:
            dwdp[8] = 1.0*E_P_1;
            break;
        case 9:
            dwdp[10] = 1.0*E*M1;
            break;
        case 10:
            dwdp[11] = 1.0*E*M;
            break;
        case 11:
            dwdp[12] = 1.0*E_M1;
            break;
        case 12:
            dwdp[13] = 1.0*E_M;
            break;
        case 13:
            dwdp[14] = 1.0*E*P;
            break;
        case 14:
            dwdp[15] = 1.0*E_P_2;
            break;
        case 15:
            dwdp[2] = 1.0*E*P21;
            break;
        case 16:
            dwdp[3] = 1.0*E*P2;
            break;
        case 17:
            dwdp[4] = 1.0*E_P21;
            break;
        case 18:
            dwdp[5] = 1.0*E_P2;
            break;
        case 19:
            dwdp[9] = 1.0*M;
            break;
        case 20:
            dwdp[6] = 1.0*E_P_1;
            break;
        case 21:
            dwdp[7] = 1.0*E_P_2;
            break;
    }
}