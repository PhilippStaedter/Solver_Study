#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model1_liu2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = C1;
    y[1] = C2;
    y[2] = C2a;
    y[3] = C2b;
    y[4] = C3;
    y[5] = C3a;
    y[6] = C3b;
    y[7] = C4;
    y[8] = C4BP;
    y[9] = C4BP_C4b;
    y[10] = C4BP_GlcNac_LF_CRP;
    y[11] = C4BP_PC_CRP;
    y[12] = C4BP_PC_CRP_LF;
    y[13] = C4a;
    y[14] = C4b;
    y[15] = C4b_C2a;
    y[16] = C4b_C2a_C4BP;
    y[17] = CRP;
    y[18] = GlcNac;
    y[19] = GlcNac_HF;
    y[20] = GlcNac_HF_MASP;
    y[21] = GlcNac_LF;
    y[22] = GlcNac_LF_C1_MASP;
    y[23] = GlcNac_LF_CRP;
    y[24] = GlcNac_LF_CRP_C1;
    y[25] = GlcNac_LF_CRP_MASP;
    y[26] = GlcNac_LF_MASP;
    y[27] = HF;
    y[28] = LF;
    y[29] = MASP;
    y[30] = PC;
    y[31] = PC_CRP;
    y[32] = PC_CRP_C1;
    y[33] = PC_CRP_LF;
    y[34] = PC_CRP_LF_C1;
    y[35] = PC_CRP_LF_C1_MASP;
    y[36] = PC_CRP_LF_MASP;
    y[37] = X;
    y[38] = dC3b;
    y[39] = dC4b_C2a;
    y[40] = dC4b_C2a_C4BP;
    y[41] = iC4b_C2a;
}