#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kirschner(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = -V*s2/pow(V + b1, 2);
            break;
        case 1:
            dwdp[4] = -V*g/pow(V + b2, 2);
            break;
        case 2:
            dwdp[5] = T*V;
            break;
        case 3:
            dwdp[4] = V/(V + b2);
            break;
        case 4:
            dwdp[3] = T*V;
            break;
        case 5:
            dwdp[2] = T;
            break;
        case 6:
            dwdp[0] = 1;
            break;
        case 7:
            dwdp[1] = V/(V + b1);
            break;
    }
}