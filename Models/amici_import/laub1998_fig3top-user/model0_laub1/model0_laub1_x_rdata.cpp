#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_laub1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = species_0;
    x_rdata[1] = species_1;
    x_rdata[2] = species_2;
    x_rdata[3] = species_3;
    x_rdata[4] = species_4;
    x_rdata[5] = species_5;
    x_rdata[6] = species_6;
}