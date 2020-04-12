#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_laub1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*species_6;
            break;
        case 1:
            dwdp[1] = 1.0*species_4;
            break;
        case 2:
            dwdp[2] = 1.0*species_4;
            break;
        case 3:
            dwdp[3] = 1.0*species_0;
            break;
        case 4:
            dwdp[4] = 1.0*species_0;
            break;
        case 5:
            dwdp[5] = 1.0*species_2*species_5;
            break;
        case 6:
            dwdp[6] = 1.0*species_1;
            break;
        case 7:
            dwdp[7] = 1.0*species_2;
            break;
        case 8:
            dwdp[8] = 1.0*species_5;
            break;
        case 9:
            dwdp[9] = 1.0*species_2*species_6;
            break;
        case 10:
            dwdp[10] = 1.0;
            break;
        case 11:
            dwdp[11] = 1.0*species_3*species_6;
            break;
        case 12:
            dwdp[12] = 1.0*species_4;
            break;
        case 13:
            dwdp[13] = 1.0*species_1*species_3;
            break;
    }
}