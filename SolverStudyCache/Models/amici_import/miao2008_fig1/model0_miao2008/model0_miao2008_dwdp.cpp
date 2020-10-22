#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_miao2008(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = T*Tm;
            break;
        case 1:
            dwdp[3] = T*Tmw;
            break;
        case 2:
            dwdp[2] = T*Tw;
            break;
        case 3:
            dwdp[5] = Tm*Tw;
            break;
        case 4:
            dwdp[7] = Tm*Tw;
            break;
        case 5:
            dwdp[0] = T;
            break;
        case 6:
            dwdp[4] = Tm;
            break;
        case 7:
            dwdp[8] = Tmw;
            break;
        case 8:
            dwdp[6] = Tw;
            break;
    }
}