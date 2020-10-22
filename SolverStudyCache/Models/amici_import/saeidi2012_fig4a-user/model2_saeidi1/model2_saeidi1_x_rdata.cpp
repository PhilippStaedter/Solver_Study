#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model2_saeidi1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = s1;
    x_rdata[1] = s16;
    x_rdata[2] = s17;
    x_rdata[3] = s19;
    x_rdata[4] = s2;
    x_rdata[5] = s3;
    x_rdata[6] = s4;
    x_rdata[7] = s42;
    x_rdata[8] = s45;
    x_rdata[9] = s5;
}