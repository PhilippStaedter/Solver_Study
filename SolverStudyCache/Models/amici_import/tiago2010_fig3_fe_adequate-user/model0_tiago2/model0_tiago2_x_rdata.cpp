#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_tiago2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = s1;
    x_rdata[1] = s10;
    x_rdata[2] = s11;
    x_rdata[3] = s12;
    x_rdata[4] = s13;
    x_rdata[5] = s14;
    x_rdata[6] = s15;
    x_rdata[7] = s16;
    x_rdata[8] = s17;
    x_rdata[9] = s2;
    x_rdata[10] = s3;
    x_rdata[11] = s4;
    x_rdata[12] = s5;
    x_rdata[13] = s6;
    x_rdata[14] = s7;
    x_rdata[15] = s8;
    x_rdata[16] = s9;
}