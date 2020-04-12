#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_valero(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = ATP*Vmapp1/(ATP + Kmapp1);
    w[1] = AMP*ATP*Vm2/(AMP*ATP + AMP*Km2ATP + ATP*Km2AMP + K);
    w[2] = ADP*Vmapp3/(ADP + Kmapp3);
    w[3] = Pyr*k4;
}