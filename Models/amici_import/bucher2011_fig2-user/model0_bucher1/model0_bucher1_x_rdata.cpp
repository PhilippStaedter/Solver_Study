#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_bucher1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ASL_b;
    x_rdata[1] = ASL_c;
    x_rdata[2] = ASL_m;
    x_rdata[3] = ASLoOH_b;
    x_rdata[4] = ASLoOH_c;
    x_rdata[5] = ASLoOH_m;
    x_rdata[6] = ASLpOH_b;
    x_rdata[7] = ASLpOH_c;
    x_rdata[8] = ASLpOH_m;
    x_rdata[9] = AS_b;
    x_rdata[10] = AS_c;
    x_rdata[11] = AS_m;
    x_rdata[12] = ASoOH_b;
    x_rdata[13] = ASoOH_c;
    x_rdata[14] = ASoOH_m;
    x_rdata[15] = ASpOH_b;
    x_rdata[16] = ASpOH_c;
    x_rdata[17] = ASpOH_m;
}