#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_alexander2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*A;
            break;
        case 1:
            dwdp[5] = 1.0*A;
            break;
        case 2:
            dwdp[3] = 1.0*G*v_max/(1.0*G + amici_k);
            break;
        case 3:
            dwdp[4] = 1.0*E;
            break;
        case 4:
            dwdp[2] = -1.0*G*v_max/pow(1.0*G + amici_k, 2);
            dwdp[3] = -1.0*G*f*v_max/pow(1.0*G + amici_k, 2);
            break;
        case 5:
            dwdp[7] = 1.0*A;
            break;
        case 6:
            dwdp[8] = 1.0*A;
            break;
        case 7:
            dwdp[10] = 1.0*E;
            break;
        case 8:
            dwdp[11] = 1.0*G;
            break;
        case 9:
            dwdp[9] = 1.0*R;
            break;
        case 10:
            dwdp[6] = 1.0*A*E;
            break;
        case 11:
            dwdp[1] = 1.0*A*R;
            break;
        case 12:
            dwdp[2] = 1.0*G/(1.0*G + amici_k);
            dwdp[3] = 1.0*G*f/(1.0*G + amici_k);
            break;
    }
}