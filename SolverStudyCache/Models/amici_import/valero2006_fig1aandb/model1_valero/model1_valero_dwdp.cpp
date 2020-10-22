#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_valero(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = -AMP*ATP*Vm2/pow(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K, 2);
            break;
        case 1:
            dwdp[1] = -AMP*pow(ATP, 2)*Vm2/pow(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K, 2);
            break;
        case 2:
            dwdp[1] = -pow(AMP, 2)*ATP*Vm2/pow(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K, 2);
            break;
        case 3:
            dwdp[0] = -ATP*Vmapp1/pow(ATP + Kmapp1, 2);
            break;
        case 4:
            dwdp[2] = -ADP*Vmapp3/pow(ADP + Kmapp3, 2);
            break;
        case 5:
            dwdp[1] = AMP*ATP/(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K);
            break;
        case 6:
            dwdp[0] = ATP/(ATP + Kmapp1);
            break;
        case 7:
            dwdp[2] = ADP/(ADP + Kmapp3);
            break;
        case 8:
            dwdp[3] = Pyr;
            break;
    }
}