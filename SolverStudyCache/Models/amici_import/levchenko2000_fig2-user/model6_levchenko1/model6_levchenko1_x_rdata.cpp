#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model6_levchenko1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = MAPK;
    x_rdata[1] = MAPKMEKpp;
    x_rdata[2] = MAPKPH;
    x_rdata[3] = MAPKp;
    x_rdata[4] = MAPKpMAPKPH;
    x_rdata[5] = MAPKpMEKpp;
    x_rdata[6] = MAPKpp;
    x_rdata[7] = MAPKppMAPKPH;
    x_rdata[8] = MEK;
    x_rdata[9] = MEKPH;
    x_rdata[10] = MEKRAFp;
    x_rdata[11] = MEKp;
    x_rdata[12] = MEKpMEKPH;
    x_rdata[13] = MEKpRAFp;
    x_rdata[14] = MEKpp;
    x_rdata[15] = MEKppMEKPH;
    x_rdata[16] = RAF;
    x_rdata[17] = RAFK;
    x_rdata[18] = RAFPH;
    x_rdata[19] = RAFRAFK;
    x_rdata[20] = RAFp;
    x_rdata[21] = RAFpRAFPH;
}