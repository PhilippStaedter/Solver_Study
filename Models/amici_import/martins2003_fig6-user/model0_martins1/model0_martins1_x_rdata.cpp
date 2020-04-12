#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_martins1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = AA;
    x_rdata[1] = Cn;
    x_rdata[2] = DFG;
    x_rdata[3] = E1;
    x_rdata[4] = E2;
    x_rdata[5] = FA;
    x_rdata[6] = Fru;
    x_rdata[7] = Glu;
    x_rdata[8] = Gly;
    x_rdata[9] = MG;
    x_rdata[10] = Man;
    x_rdata[11] = Mel;
    x_rdata[12] = _1DG;
    x_rdata[13] = _3DG;
}