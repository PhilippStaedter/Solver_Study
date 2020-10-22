#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model1_levchenko2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = C1;
    x_rdata[1] = C2;
    x_rdata[2] = C3;
    x_rdata[3] = C4;
    x_rdata[4] = C5;
    x_rdata[5] = C6;
    x_rdata[6] = C7;
    x_rdata[7] = C8;
    x_rdata[8] = C9;
    x_rdata[9] = MAPK;
    x_rdata[10] = MAPKMEKpp;
    x_rdata[11] = MAPKPH;
    x_rdata[12] = MAPKp;
    x_rdata[13] = MAPKpMAPKPH;
    x_rdata[14] = MAPKpMEKpp;
    x_rdata[15] = MAPKpp;
    x_rdata[16] = MAPKppMAPKPH;
    x_rdata[17] = MEK;
    x_rdata[18] = MEKPH;
    x_rdata[19] = MEKRAFp;
    x_rdata[20] = MEKp;
    x_rdata[21] = MEKpMEKPH;
    x_rdata[22] = MEKpRAFp;
    x_rdata[23] = MEKpp;
    x_rdata[24] = MEKppMEKPH;
    x_rdata[25] = RAF;
    x_rdata[26] = RAFK;
    x_rdata[27] = RAFPH;
    x_rdata[28] = RAFRAFK;
    x_rdata[29] = RAFp;
    x_rdata[30] = RAFpRAFPH;
}