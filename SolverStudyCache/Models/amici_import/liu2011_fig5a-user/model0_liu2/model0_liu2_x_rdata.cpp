#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_liu2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = C1;
    x_rdata[1] = C2;
    x_rdata[2] = C2a;
    x_rdata[3] = C2b;
    x_rdata[4] = C3;
    x_rdata[5] = C3a;
    x_rdata[6] = C3b;
    x_rdata[7] = C4;
    x_rdata[8] = C4BP;
    x_rdata[9] = C4BP_C4b;
    x_rdata[10] = C4BP_GlcNac_LF_CRP;
    x_rdata[11] = C4BP_PC_CRP;
    x_rdata[12] = C4BP_PC_CRP_LF;
    x_rdata[13] = C4a;
    x_rdata[14] = C4b;
    x_rdata[15] = C4b_C2a;
    x_rdata[16] = C4b_C2a_C4BP;
    x_rdata[17] = CRP;
    x_rdata[18] = GlcNac;
    x_rdata[19] = GlcNac_HF;
    x_rdata[20] = GlcNac_HF_MASP;
    x_rdata[21] = GlcNac_LF;
    x_rdata[22] = GlcNac_LF_C1_MASP;
    x_rdata[23] = GlcNac_LF_CRP;
    x_rdata[24] = GlcNac_LF_CRP_C1;
    x_rdata[25] = GlcNac_LF_CRP_MASP;
    x_rdata[26] = GlcNac_LF_MASP;
    x_rdata[27] = HF;
    x_rdata[28] = LF;
    x_rdata[29] = MASP;
    x_rdata[30] = PC;
    x_rdata[31] = PC_CRP;
    x_rdata[32] = PC_CRP_C1;
    x_rdata[33] = PC_CRP_LF;
    x_rdata[34] = PC_CRP_LF_C1;
    x_rdata[35] = PC_CRP_LF_C1_MASP;
    x_rdata[36] = PC_CRP_LF_MASP;
    x_rdata[37] = X;
    x_rdata[38] = dC3b;
    x_rdata[39] = dC4b_C2a;
    x_rdata[40] = dC4b_C2a_C4BP;
    x_rdata[41] = iC4b_C2a;
}