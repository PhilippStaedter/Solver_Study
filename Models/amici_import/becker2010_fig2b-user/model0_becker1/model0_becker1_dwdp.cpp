#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_becker1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*kt;
            break;
        case 1:
            dwdp[7] = 1.0*Epo_EpoRi;
            break;
        case 2:
            dwdp[6] = 1.0*Epo_EpoRi;
            break;
        case 3:
            dwdp[4] = 1.0*Epo_EpoR;
            break;
        case 4:
            dwdp[5] = 1.0*Epo_EpoRi;
            break;
        case 5:
            dwdp[3] = 1.0*Epo_EpoR;
            break;
        case 6:
            dwdp[2] = 1.0*Epo*EpoR;
            break;
        case 7:
            dwdp[0] = 1.0*Bmax;
            dwdp[1] = 1.0*EpoR;
            break;
    }
}