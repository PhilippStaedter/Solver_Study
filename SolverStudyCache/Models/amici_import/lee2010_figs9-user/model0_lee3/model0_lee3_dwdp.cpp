#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_lee3(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*II;
            break;
        case 1:
            dwdp[1] = 1.0*M;
            break;
        case 2:
            dwdp[2] = 1.0*II;
            break;
        case 3:
            dwdp[3] = 1.0*P2;
            break;
    }
}