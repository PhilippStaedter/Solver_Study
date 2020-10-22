#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_esteban2013(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[10] = -M*V*rhoM/pow(CM + V, 2);
            break;
        case 1:
            dwdp[5] = -T*V*rhoT/pow(CT + V, 2);
            break;
        case 2:
            dwdp[12] = M;
            break;
        case 3:
            dwdp[1] = MStar;
            break;
        case 4:
            dwdp[7] = T;
            break;
        case 5:
            dwdp[8] = TStar;
            break;
        case 6:
            dwdp[4] = V;
            break;
        case 7:
            dwdp[11] = M*V;
            break;
        case 8:
            dwdp[6] = T*V;
            break;
        case 9:
            dwdp[3] = MStar;
            break;
        case 10:
            dwdp[2] = TStar;
            break;
        case 11:
            dwdp[10] = M*V/(CM + V);
            break;
        case 12:
            dwdp[5] = T*V/(CT + V);
            break;
        case 13:
            dwdp[9] = 1;
            break;
        case 14:
            dwdp[0] = 1;
            break;
    }
}