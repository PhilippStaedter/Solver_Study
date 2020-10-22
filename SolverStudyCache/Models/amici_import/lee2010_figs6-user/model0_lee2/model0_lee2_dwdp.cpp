#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_lee2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*E_P_1;
            break;
        case 1:
            dwdp[4] = -1.0*E_M;
            break;
        case 2:
            dwdp[6] = -1.0*E_P_2;
            break;
        case 3:
            dwdp[1] = -1.0*E_P2;
            break;
        case 4:
            dwdp[0] = 1.0*E*P;
            break;
        case 5:
            dwdp[3] = 1.0*E_P_1;
            break;
        case 6:
            dwdp[4] = 1.0*E*M;
            break;
        case 7:
            dwdp[5] = 1.0*E_M;
            break;
        case 8:
            dwdp[6] = 1.0*E*P;
            break;
        case 9:
            dwdp[7] = 1.0*E_P_2;
            break;
        case 10:
            dwdp[1] = 1.0*E*P2;
            break;
        case 11:
            dwdp[2] = 1.0*E_P2;
            break;
    }
}