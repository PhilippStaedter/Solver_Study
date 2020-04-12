#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model1_hald(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ACA;
    x_rdata[1] = ACA0;
    x_rdata[2] = ADP;
    x_rdata[3] = AMP;
    x_rdata[4] = ATP;
    x_rdata[5] = DHAP;
    x_rdata[6] = DHAPCN;
    x_rdata[7] = DPG;
    x_rdata[8] = EtOH;
    x_rdata[9] = EtOH0;
    x_rdata[10] = F6P;
    x_rdata[11] = FBP;
    x_rdata[12] = G6P;
    x_rdata[13] = GAP;
    x_rdata[14] = Glc;
    x_rdata[15] = Glc0;
    x_rdata[16] = Glyc;
    x_rdata[17] = Glyc0;
    x_rdata[18] = HCN;
    x_rdata[19] = HCN0;
    x_rdata[20] = NAD;
    x_rdata[21] = NADH;
    x_rdata[22] = OAc;
    x_rdata[23] = OAc0;
    x_rdata[24] = PEP;
    x_rdata[25] = Pyr;
    x_rdata[26] = PyrCN;
    x_rdata[27] = X;
    x_rdata[28] = drain;
    x_rdata[29] = glycogen;
    x_rdata[30] = lacto;
}