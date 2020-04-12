#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Liu2011(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CRP;
    x_rdata[1] = PC;
    x_rdata[2] = PC_CRP;
    x_rdata[3] = C4;
    x_rdata[4] = C4a;
    x_rdata[5] = C4b;
    x_rdata[6] = C2;
    x_rdata[7] = C1;
    x_rdata[8] = PC_CRP_C1;
    x_rdata[9] = C2a;
    x_rdata[10] = C2b;
    x_rdata[11] = C4b_C2a;
    x_rdata[12] = C3;
    x_rdata[13] = C3a;
    x_rdata[14] = C3b;
    x_rdata[15] = dC3b;
    x_rdata[16] = MASP;
    x_rdata[17] = dC4b_C2a;
    x_rdata[18] = GlcNac;
    x_rdata[19] = GlcNac_LF;
    x_rdata[20] = LF;
    x_rdata[21] = GlcNac_LF_MASP;
    x_rdata[22] = PC_CRP_LF;
    x_rdata[23] = PC_CRP_LF_MASP;
    x_rdata[24] = GlcNac_LF_CRP;
    x_rdata[25] = GlcNac_LF_CRP_C1;
    x_rdata[26] = C4BP;
    x_rdata[27] = C4BP_PC_CRP;
    x_rdata[28] = C4BP_GlcNac_LF_CRP;
    x_rdata[29] = iC4b_C2a;
    x_rdata[30] = C4BP_C4b;
    x_rdata[31] = C4b_C2a_C4BP;
    x_rdata[32] = dC4b_C2a_C4BP;
    x_rdata[33] = PC_CRP_LF_C1;
    x_rdata[34] = C4BP_PC_CRP_LF;
    x_rdata[35] = GlcNac_LF_CRP_MASP;
    x_rdata[36] = PC_CRP_LF_C1_MASP;
    x_rdata[37] = GlcNac_LF_C1_MASP;
    x_rdata[38] = GlcNac_HF;
    x_rdata[39] = HF;
    x_rdata[40] = GlcNac_HF_MASP;
    x_rdata[41] = X;
}