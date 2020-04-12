#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kholodenko2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = MAPK;
    x_rdata[1] = MAPK_P;
    x_rdata[2] = MAPK_PP;
    x_rdata[3] = MKK;
    x_rdata[4] = MKKK;
    x_rdata[5] = MKKK_P;
    x_rdata[6] = MKK_P;
    x_rdata[7] = MKK_PP;
}