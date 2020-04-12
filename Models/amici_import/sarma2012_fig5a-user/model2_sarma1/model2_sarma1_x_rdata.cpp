#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model2_sarma1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = species_1;
    x_rdata[1] = species_10;
    x_rdata[2] = species_11;
    x_rdata[3] = species_12;
    x_rdata[4] = species_13;
    x_rdata[5] = species_14;
    x_rdata[6] = species_15;
    x_rdata[7] = species_16;
    x_rdata[8] = species_17;
    x_rdata[9] = species_18;
    x_rdata[10] = species_19;
    x_rdata[11] = species_2;
    x_rdata[12] = species_20;
    x_rdata[13] = species_21;
    x_rdata[14] = species_22;
    x_rdata[15] = species_23;
    x_rdata[16] = species_24;
    x_rdata[17] = species_25;
    x_rdata[18] = species_26;
    x_rdata[19] = species_27;
    x_rdata[20] = species_3;
    x_rdata[21] = species_4;
    x_rdata[22] = species_5;
    x_rdata[23] = species_6;
    x_rdata[24] = species_7;
    x_rdata[25] = species_8;
    x_rdata[26] = species_9;
}