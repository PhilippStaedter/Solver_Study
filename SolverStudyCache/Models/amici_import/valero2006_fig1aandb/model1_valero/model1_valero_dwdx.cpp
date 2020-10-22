#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model1_valero(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -ADP*Vmapp3/pow(ADP + Kmapp3, 2) + Vmapp3/(ADP + Kmapp3);
    dwdx[1] = AMP*ATP*Vm2*(-ATP - Km2ATP)/pow(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K, 2) + ATP*Vm2/(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K);
    dwdx[2] = -ATP*Vmapp1/pow(ATP + Kmapp1, 2) + Vmapp1/(ATP + Kmapp1);
    dwdx[3] = AMP*ATP*Vm2*(-AMP - Km2AMP)/pow(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K, 2) + AMP*Vm2/(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K);
    dwdx[4] = k4;
}