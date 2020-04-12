#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Yang2007(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = x1;
    x_rdata[1] = x10;
    x_rdata[2] = x11;
    x_rdata[3] = x12;
    x_rdata[4] = x13;
    x_rdata[5] = x14;
    x_rdata[6] = x15;
    x_rdata[7] = x16;
    x_rdata[8] = x17;
    x_rdata[9] = x18;
    x_rdata[10] = x19;
    x_rdata[11] = x2;
    x_rdata[12] = x20;
    x_rdata[13] = x21;
    x_rdata[14] = x22;
    x_rdata[15] = x23;
    x_rdata[16] = x24;
    x_rdata[17] = x25;
    x_rdata[18] = x3;
    x_rdata[19] = x4;
    x_rdata[20] = x5;
    x_rdata[21] = x6;
    x_rdata[22] = x7;
    x_rdata[23] = x8;
    x_rdata[24] = x9;
}