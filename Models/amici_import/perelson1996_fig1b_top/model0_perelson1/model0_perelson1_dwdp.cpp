#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_perelson1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[4] = Tstar*delta;
            break;
        case 1:
            dwdp[0] = Vin*amici_k;
            break;
        case 2:
            dwdp[2] = Vin;
            dwdp[3] = Vni;
            break;
        case 3:
            dwdp[1] = Tstar;
            dwdp[4] = NN*Tstar;
            break;
        case 4:
            dwdp[0] = T0*Vin;
            break;
    }
}