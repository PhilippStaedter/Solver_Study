#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model1_mellor1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = species_1;
    x_rdata[1] = species_10;
    x_rdata[2] = species_11;
    x_rdata[3] = species_12;
    x_rdata[4] = species_13;
    x_rdata[5] = species_14;
    x_rdata[6] = species_15;
    x_rdata[7] = species_7;
    x_rdata[8] = species_8;
    x_rdata[9] = species_9;
}