#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_ma(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = CAR1;
            break;
        case 1:
            dwdp[1] = REGA*incAMP;
            break;
        case 2:
            dwdp[2] = ACA;
            break;
        case 3:
            dwdp[3] = excAMP;
            break;
        case 4:
            dwdp[4] = excAMP;
            break;
        case 5:
            dwdp[5] = CAR1;
            break;
        case 6:
            dwdp[6] = ACA*PKA;
            break;
        case 7:
            dwdp[7] = incAMP;
            break;
        case 8:
            dwdp[8] = PKA;
            break;
        case 9:
            dwdp[9] = CAR1;
            break;
        case 10:
            dwdp[10] = ERK2*PKA;
            break;
        case 11:
            dwdp[11] = 1;
            break;
        case 12:
            dwdp[12] = ERK2*REGA;
            break;
        case 13:
            dwdp[13] = ACA;
            break;
    }
}