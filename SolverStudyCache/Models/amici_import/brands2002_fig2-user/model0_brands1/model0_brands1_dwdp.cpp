#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_brands1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = Glu;
            break;
        case 1:
            dwdp[1] = Fru*lys_R;
            break;
        case 2:
            dwdp[2] = AMP;
            break;
        case 3:
            dwdp[3] = Fru;
            break;
        case 4:
            dwdp[4] = Glu;
            break;
        case 5:
            dwdp[5] = Fru;
            break;
        case 6:
            dwdp[6] = Fru;
            break;
        case 7:
            dwdp[7] = Triose;
            break;
        case 8:
            dwdp[8] = Glu*lys_R;
            break;
        case 9:
            dwdp[9] = Amadori;
            break;
        case 10:
            dwdp[10] = Amadori;
            break;
    }
}