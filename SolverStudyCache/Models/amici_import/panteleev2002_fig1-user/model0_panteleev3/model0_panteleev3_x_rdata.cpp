#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_panteleev3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = TFPI;
    x_rdata[1] = VIIa_TF;
    x_rdata[2] = VIIa_TF_X;
    x_rdata[3] = VIIa_TF_Xa;
    x_rdata[4] = X;
    x_rdata[5] = Xa;
    x_rdata[6] = Xa_TFPI;
    x_rdata[7] = Xa_TFPI_VIIa_TF;
}