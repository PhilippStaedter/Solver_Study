#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_brands1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = AMP;
    x_rdata[1] = Acetic_acid;
    x_rdata[2] = Amadori;
    x_rdata[3] = C5;
    x_rdata[4] = Cn;
    x_rdata[5] = Formic_acid;
    x_rdata[6] = Fru;
    x_rdata[7] = Glu;
    x_rdata[8] = Melanoidin;
    x_rdata[9] = Triose;
    x_rdata[10] = lys_R;
}