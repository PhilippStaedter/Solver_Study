#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Liu2011(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CRP;
    y[1] = PC;
    y[2] = PC_CRP;
    y[3] = C4;
    y[4] = C4a;
    y[5] = C4b;
    y[6] = C2;
    y[7] = C1;
    y[8] = PC_CRP_C1;
    y[9] = C2a;
    y[10] = C2b;
    y[11] = C4b_C2a;
    y[12] = C3;
    y[13] = C3a;
    y[14] = C3b;
    y[15] = dC3b;
    y[16] = MASP;
    y[17] = dC4b_C2a;
    y[18] = GlcNac;
    y[19] = GlcNac_LF;
    y[20] = LF;
    y[21] = GlcNac_LF_MASP;
    y[22] = PC_CRP_LF;
    y[23] = PC_CRP_LF_MASP;
    y[24] = GlcNac_LF_CRP;
    y[25] = GlcNac_LF_CRP_C1;
    y[26] = C4BP;
    y[27] = C4BP_PC_CRP;
    y[28] = C4BP_GlcNac_LF_CRP;
    y[29] = iC4b_C2a;
    y[30] = C4BP_C4b;
    y[31] = C4b_C2a_C4BP;
    y[32] = dC4b_C2a_C4BP;
    y[33] = PC_CRP_LF_C1;
    y[34] = C4BP_PC_CRP_LF;
    y[35] = GlcNac_LF_CRP_MASP;
    y[36] = PC_CRP_LF_C1_MASP;
    y[37] = GlcNac_LF_C1_MASP;
    y[38] = GlcNac_HF;
    y[39] = HF;
    y[40] = GlcNac_HF_MASP;
    y[41] = X;
}