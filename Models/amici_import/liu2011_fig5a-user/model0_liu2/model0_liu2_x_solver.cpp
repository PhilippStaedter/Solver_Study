#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_liu2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = C1;
    x_solver[1] = C2;
    x_solver[2] = C2a;
    x_solver[3] = C2b;
    x_solver[4] = C3;
    x_solver[5] = C3a;
    x_solver[6] = C3b;
    x_solver[7] = C4;
    x_solver[8] = C4BP;
    x_solver[9] = C4BP_C4b;
    x_solver[10] = C4BP_GlcNac_LF_CRP;
    x_solver[11] = C4BP_PC_CRP;
    x_solver[12] = C4BP_PC_CRP_LF;
    x_solver[13] = C4a;
    x_solver[14] = C4b;
    x_solver[15] = C4b_C2a;
    x_solver[16] = C4b_C2a_C4BP;
    x_solver[17] = CRP;
    x_solver[18] = GlcNac;
    x_solver[19] = GlcNac_HF;
    x_solver[20] = GlcNac_HF_MASP;
    x_solver[21] = GlcNac_LF;
    x_solver[22] = GlcNac_LF_C1_MASP;
    x_solver[23] = GlcNac_LF_CRP;
    x_solver[24] = GlcNac_LF_CRP_C1;
    x_solver[25] = GlcNac_LF_CRP_MASP;
    x_solver[26] = GlcNac_LF_MASP;
    x_solver[27] = HF;
    x_solver[28] = LF;
    x_solver[29] = MASP;
    x_solver[30] = PC;
    x_solver[31] = PC_CRP;
    x_solver[32] = PC_CRP_C1;
    x_solver[33] = PC_CRP_LF;
    x_solver[34] = PC_CRP_LF_C1;
    x_solver[35] = PC_CRP_LF_C1_MASP;
    x_solver[36] = PC_CRP_LF_MASP;
    x_solver[37] = X;
    x_solver[38] = dC3b;
    x_solver[39] = dC4b_C2a;
    x_solver[40] = dC4b_C2a_C4BP;
    x_solver[41] = iC4b_C2a;
}