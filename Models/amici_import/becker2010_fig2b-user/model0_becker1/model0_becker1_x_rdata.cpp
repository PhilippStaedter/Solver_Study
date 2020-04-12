#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_becker1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Epo;
    x_rdata[1] = EpoR;
    x_rdata[2] = Epo_EpoR;
    x_rdata[3] = Epo_EpoRi;
    x_rdata[4] = dEpoe;
    x_rdata[5] = dEpoi;
}