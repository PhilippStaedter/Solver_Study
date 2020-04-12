#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_fung1_Fig3A_Vgly_0_5(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = AcCoA;
    x_rdata[1] = AcP;
    x_rdata[2] = Acs;
    x_rdata[3] = HOAc;
    x_rdata[4] = HOAc_E;
    x_rdata[5] = LacI;
    x_rdata[6] = OAc;
    x_rdata[7] = Pta;
}