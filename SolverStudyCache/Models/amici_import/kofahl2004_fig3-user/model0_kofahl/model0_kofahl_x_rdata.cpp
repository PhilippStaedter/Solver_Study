#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_kofahl(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Bar1;
    x_rdata[1] = Bar1a;
    x_rdata[2] = Bar1aex;
    x_rdata[3] = Cdc28;
    x_rdata[4] = Far1;
    x_rdata[5] = Far1PP;
    x_rdata[6] = Far1U;
    x_rdata[7] = Fus3;
    x_rdata[8] = Fus3PP;
    x_rdata[9] = GaGDP;
    x_rdata[10] = GaGTP;
    x_rdata[11] = Gabc;
    x_rdata[12] = Gbc;
    x_rdata[13] = Sst2;
    x_rdata[14] = Ste11;
    x_rdata[15] = Ste12;
    x_rdata[16] = Ste12a;
    x_rdata[17] = Ste2;
    x_rdata[18] = Ste20;
    x_rdata[19] = Ste2a;
    x_rdata[20] = Ste5;
    x_rdata[21] = Ste7;
    x_rdata[22] = alpha;
    x_rdata[23] = complexA;
    x_rdata[24] = complexB;
    x_rdata[25] = complexC;
    x_rdata[26] = complexD;
    x_rdata[27] = complexE;
    x_rdata[28] = complexF;
    x_rdata[29] = complexG;
    x_rdata[30] = complexH;
    x_rdata[31] = complexI;
    x_rdata[32] = complexK;
    x_rdata[33] = complexL;
    x_rdata[34] = complexM;
    x_rdata[35] = complexN;
    x_rdata[36] = amici_p;
}