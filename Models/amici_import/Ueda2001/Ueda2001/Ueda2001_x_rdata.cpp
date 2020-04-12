#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Ueda2001(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EmptySet;
    x_rdata[1] = CCc;
    x_rdata[2] = CCn;
    x_rdata[3] = Clkc;
    x_rdata[4] = Clkm;
    x_rdata[5] = Perc;
    x_rdata[6] = Perm;
    x_rdata[7] = PTc;
    x_rdata[8] = PTn;
    x_rdata[9] = Timc;
    x_rdata[10] = Timm;
    x_rdata[11] = species_0000012;
    x_rdata[12] = species_0000013;
}