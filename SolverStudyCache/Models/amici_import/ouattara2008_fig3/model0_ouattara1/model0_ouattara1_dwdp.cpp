#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_ouattara1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[2] = T*V;
            break;
        case 1:
            dwdp[5] = V;
            break;
        case 2:
            dwdp[1] = T;
            break;
        case 3:
            dwdp[4] = TStar;
            break;
        case 4:
            dwdp[3] = TStar;
            break;
        case 5:
            dwdp[0] = 1;
            break;
    }
}