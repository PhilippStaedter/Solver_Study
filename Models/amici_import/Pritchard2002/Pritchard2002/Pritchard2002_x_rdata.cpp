#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Pritchard2002(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = GLCo;
    x_rdata[1] = GLCi;
    x_rdata[2] = ATP;
    x_rdata[3] = G6P;
    x_rdata[4] = ADP;
    x_rdata[5] = F6P;
    x_rdata[6] = F16bP;
    x_rdata[7] = AMP;
    x_rdata[8] = F26bP;
    x_rdata[9] = DHAP;
    x_rdata[10] = GAP;
    x_rdata[11] = NAD;
    x_rdata[12] = BPG;
    x_rdata[13] = NADH;
    x_rdata[14] = P3G;
    x_rdata[15] = P2G;
    x_rdata[16] = PEP;
    x_rdata[17] = PYR;
    x_rdata[18] = AcAld;
    x_rdata[19] = CO2;
    x_rdata[20] = EtOH;
    x_rdata[21] = Glycerol;
    x_rdata[22] = Glycogen;
    x_rdata[23] = Trehalose;
    x_rdata[24] = Succinate;
}