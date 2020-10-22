#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model2_(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EtOH;
    x_rdata[1] = Glycerol;
    x_rdata[2] = PiVac;
    x_rdata[3] = atp;
    x_rdata[4] = fbp;
    x_rdata[5] = g6p;
    x_rdata[6] = phos;
}