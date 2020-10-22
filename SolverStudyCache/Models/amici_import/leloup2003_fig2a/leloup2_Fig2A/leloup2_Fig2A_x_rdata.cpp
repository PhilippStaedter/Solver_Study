#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_leloup2_Fig2A(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = species_0;
    x_rdata[1] = species_1;
    x_rdata[2] = species_10;
    x_rdata[3] = species_11;
    x_rdata[4] = species_12;
    x_rdata[5] = species_13;
    x_rdata[6] = species_14;
    x_rdata[7] = species_15;
    x_rdata[8] = species_2;
    x_rdata[9] = species_3;
    x_rdata[10] = species_4;
    x_rdata[11] = species_5;
    x_rdata[12] = species_6;
    x_rdata[13] = species_7;
    x_rdata[14] = species_8;
    x_rdata[15] = species_9;
}