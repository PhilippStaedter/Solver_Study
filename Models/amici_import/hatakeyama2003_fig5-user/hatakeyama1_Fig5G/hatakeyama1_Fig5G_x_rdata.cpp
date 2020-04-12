#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_hatakeyama1_Fig5G(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Akt;
    x_rdata[1] = AktPIP;
    x_rdata[2] = AktPIP3;
    x_rdata[3] = AktPIPP;
    x_rdata[4] = E;
    x_rdata[5] = ERK;
    x_rdata[6] = ERKP;
    x_rdata[7] = ERKPP;
    x_rdata[8] = GS;
    x_rdata[9] = HRG;
    x_rdata[10] = MEK;
    x_rdata[11] = MEKP;
    x_rdata[12] = MEKPP;
    x_rdata[13] = MKP3;
    x_rdata[14] = PI3K;
    x_rdata[15] = PI3Kstar;
    x_rdata[16] = PIP3;
    x_rdata[17] = PP2A;
    x_rdata[18] = P_I;
    x_rdata[19] = R;
    x_rdata[20] = RHRG;
    x_rdata[21] = RHRG2;
    x_rdata[22] = RP;
    x_rdata[23] = RPI3K;
    x_rdata[24] = RPI3Kstar;
    x_rdata[25] = RShGS;
    x_rdata[26] = RShP;
    x_rdata[27] = RShc;
    x_rdata[28] = Raf;
    x_rdata[29] = Rafstar;
    x_rdata[30] = RasGDP;
    x_rdata[31] = RasGTP;
    x_rdata[32] = ShGS;
    x_rdata[33] = ShP;
    x_rdata[34] = Shc;
    x_rdata[35] = internalization;
}