#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_piedrafita1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*S*STU;
            break;
        case 1:
            dwdp[1] = 1.0*STUSU;
            break;
        case 2:
            dwdp[1] = -1.0*STU*SU;
            break;
        case 3:
            dwdp[0] = -1.0*STUS;
            break;
        case 4:
            dwdp[3] = 1.0*STUS*T;
            break;
        case 5:
            dwdp[3] = -1.0*STUST;
            break;
        case 6:
            dwdp[4] = 1.0*STUST;
            break;
        case 7:
            dwdp[4] = -1.0*ST*STU;
            break;
        case 8:
            dwdp[2] = 1.0*ST;
            dwdp[5] = 1.0*STU;
            dwdp[9] = 1.0*SU;
            break;
        case 9:
            dwdp[6] = 1.0*ST*SU;
            break;
        case 10:
            dwdp[6] = -1.0*SUST;
            break;
        case 11:
            dwdp[7] = 1.0*SUST*U;
            break;
        case 12:
            dwdp[7] = -1.0*SUSTU;
            break;
        case 13:
            dwdp[8] = 1.0*SUSTU;
            break;
        case 14:
            dwdp[8] = -1.0*STU*SU;
            break;
        case 15:
            dwdp[10] = 1.0*STUS*U;
            break;
        case 16:
            dwdp[10] = -1.0*STUSU;
            break;
    }
}