#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_martins1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = DFG;
            break;
        case 1:
            dwdp[1] = E1;
            break;
        case 2:
            dwdp[2] = E1;
            break;
        case 3:
            dwdp[3] = Man;
            break;
        case 4:
            dwdp[4] = Glu;
            break;
        case 5:
            dwdp[5] = Cn*Gly;
            break;
        case 6:
            dwdp[6] = Cn;
            break;
        case 7:
            dwdp[7] = E2;
            break;
        case 8:
            dwdp[8] = DFG;
            break;
        case 9:
            dwdp[9] = DFG;
            break;
        case 10:
            dwdp[10] = E1;
            break;
        case 11:
            dwdp[11] = _3DG;
            break;
        case 12:
            dwdp[12] = _3DG;
            break;
        case 13:
            dwdp[13] = E2;
            break;
        case 14:
            dwdp[14] = _1DG;
            break;
        case 15:
            dwdp[15] = _1DG;
            break;
    }
}