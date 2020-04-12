#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_ferreira1_Fig2F(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Amadori;
    x_rdata[1] = CML;
    x_rdata[2] = Glucose;
    x_rdata[3] = Glyoxal;
    x_rdata[4] = Lysine;
    x_rdata[5] = Schiff;
}