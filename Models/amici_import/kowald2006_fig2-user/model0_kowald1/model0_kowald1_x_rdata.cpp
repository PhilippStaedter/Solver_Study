#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kowald1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = species_0000001;
    x_rdata[1] = species_0000002;
    x_rdata[2] = species_0000006;
    x_rdata[3] = species_0000007;
    x_rdata[4] = species_0000008;
    x_rdata[5] = species_0000009;
    x_rdata[6] = species_0000011;
    x_rdata[7] = species_0000016;
    x_rdata[8] = species_0000017;
}