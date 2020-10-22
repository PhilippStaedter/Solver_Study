#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_zi1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = LRC_Cave;
    x_rdata[1] = LRC_EE;
    x_rdata[2] = LRC_Surf;
    x_rdata[3] = Smad2c;
    x_rdata[4] = Smad2n;
    x_rdata[5] = Smad4c;
    x_rdata[6] = Smad4n;
    x_rdata[7] = Smads_Complex_c;
    x_rdata[8] = Smads_Complex_n;
    x_rdata[9] = T1R_Cave;
    x_rdata[10] = T1R_EE;
    x_rdata[11] = T1R_Surf;
    x_rdata[12] = T2R_Cave;
    x_rdata[13] = T2R_EE;
    x_rdata[14] = T2R_Surf;
    x_rdata[15] = TGF_beta;
}