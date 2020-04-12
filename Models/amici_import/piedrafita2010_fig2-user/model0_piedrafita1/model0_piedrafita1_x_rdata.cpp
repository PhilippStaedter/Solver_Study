#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_piedrafita1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = S;
    x_rdata[1] = ST;
    x_rdata[2] = STU;
    x_rdata[3] = STUS;
    x_rdata[4] = STUST;
    x_rdata[5] = STUSU;
    x_rdata[6] = SU;
    x_rdata[7] = SUST;
    x_rdata[8] = SUSTU;
    x_rdata[9] = T;
    x_rdata[10] = U;
}