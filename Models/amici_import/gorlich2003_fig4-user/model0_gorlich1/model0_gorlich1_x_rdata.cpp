#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_gorlich1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = GDP;
    x_rdata[1] = GTP;
    x_rdata[2] = RCC1;
    x_rdata[3] = RCC1_Ran;
    x_rdata[4] = RCC1_RanGDP;
    x_rdata[5] = RCC1_RanGTP;
    x_rdata[6] = RanBP1;
    x_rdata[7] = RanGAP;
    x_rdata[8] = RanGDP_cy;
    x_rdata[9] = RanGDP_nuc;
    x_rdata[10] = RanGTP_RanBP1;
    x_rdata[11] = RanGTP_cy;
    x_rdata[12] = RanGTP_nuc;
}