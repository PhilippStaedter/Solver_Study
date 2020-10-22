#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model2_liu2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*CRP*PC*ka01_1 - 1.0*PC_CRP*ka01_2;
    w[1] = 1.0*C1*PC_CRP*ka02_1 - 1.0*PC_CRP_C1*ka02_2;
    w[2] = 1.0*C4*PC_CRP_C1*ka03_1/(C4 + ka03_2);
    w[3] = 1.0*C2*PC_CRP_C1*ka04_1/(C2 + ka04_2);
    w[4] = 1.0*GlcNac*LF*kb01_1 - 1.0*GlcNac_LF*kb01_2;
    w[5] = 1.0*GlcNac_LF*MASP*kb02_1 - 1.0*GlcNac_LF_MASP*kb02_2;
    w[6] = 1.0*C4*GlcNac_LF_MASP*kb03_1/(C4 + kb03_2);
    w[7] = 1.0*C2*GlcNac_LF_MASP*kb04_1/(C2 + kb04_2);
    w[8] = 1.0*C2a*C4b*kc01_1 - 1.0*C4b_C2a*kc01_2;
    w[9] = 1.0*C3*C4b_C2a*kc02;
    w[10] = 1.0*C3b*kc03_1 - 1.0*dC3b*kc03_2;
    w[11] = 1.0*C4b_C2a*kc04_1 - 1.0*dC4b_C2a*kc04_2;
    w[12] = 1.0*LF*PC_CRP*kd01_1 - 1.0*PC_CRP_LF*kd01_2;
    w[13] = 1.0*MASP*PC_CRP_LF*kd02_1 - 1.0*PC_CRP_LF_MASP*kd02_2;
    w[14] = 1.0*C4*PC_CRP_LF_MASP*kd03_1/(C4 + kd03_2);
    w[15] = 1.0*C2*PC_CRP_LF_MASP*kd04_1/(C2 + kd04_2);
    w[16] = 1.0*C1*PC_CRP_LF*kd05_1 - 1.0*PC_CRP_LF_C1*kd05_2;
    w[17] = 1.0*C4*PC_CRP_LF_C1*kd06_1/(C4 + kd06_2);
    w[18] = 1.0*C2*PC_CRP_LF_C1*kd07_1/(C2 + kd07_2);
    w[19] = 1.0*MASP*PC_CRP_LF_C1*kd08_1 - 1.0*PC_CRP_LF_C1_MASP*kd08_2;
    w[20] = 1.0*C1*PC_CRP_LF_MASP*kd09_1 - 1.0*PC_CRP_LF_C1_MASP*kd09_2;
    w[21] = 1.0*C4*PC_CRP_LF_C1_MASP*kd10_1/(C4 + kd10_2);
    w[22] = 1.0*C2*PC_CRP_LF_C1_MASP*kd11_1/(C2 + kd11_2);
    w[23] = 1.0*CRP*GlcNac_LF*ke01_1 - 1.0*GlcNac_LF_CRP*ke01_2;
    w[24] = 1.0*C1*GlcNac_LF_CRP*ke02_1 - 1.0*GlcNac_LF_CRP_C1*ke02_2;
    w[25] = 1.0*C4*GlcNac_LF_CRP_C1*ke03_1/(C4 + ke03_2);
    w[26] = 1.0*C2*GlcNac_LF_CRP_C1*ke04_1/(C2 + ke04_2);
    w[27] = 1.0*GlcNac_LF_CRP*MASP*ke05_1 - 1.0*GlcNac_LF_CRP_MASP*ke05_2;
    w[28] = 1.0*C4*GlcNac_LF_CRP_MASP*ke06_1/(C4 + ke06_2);
    w[29] = 1.0*C2*GlcNac_LF_CRP_MASP*ke07_1/(C2 + ke07_2);
    w[30] = 1.0*C4BP*PC_CRP*kf01_1 - 1.0*C4BP_PC_CRP*kf01_2;
    w[31] = 1.0*C4BP*GlcNac_LF_CRP*kf02_1 - 1.0*C4BP_GlcNac_LF_CRP*kf02_2;
    w[32] = 1.0*C4BP*C4b_C2a*kf03;
    w[33] = 1.0*C4BP*C4b*kf04_1 - 1.0*C4BP_C4b*kf04_2;
    w[34] = 1.0*C4BP*C4b_C2a*kf05;
    w[35] = 1.0*C4BP*C4b_C2a*kf06_1 - 1.0*C4b_C2a_C4BP*kf06_2;
    w[36] = 1.0*C4BP*dC4b_C2a*kf07_1 - 1.0*dC4b_C2a_C4BP*kf07_2;
    w[37] = -1.0*GlcNac_HF*kg01_2 + 1.0*HF*X*kg01_1;
    w[38] = 1.0*GlcNac_HF*MASP*kg02_1 - 1.0*GlcNac_HF_MASP*kg02_2;
    w[39] = 1.0*C4*GlcNac_HF_MASP*kg03_1/(C4 + kg03_2);
    w[40] = 1.0*C2*GlcNac_HF_MASP*kg04_1/(C2 + kg04_2);
    w[41] = 1.0*C4BP*t_01_k1_4;
    w[42] = 1.0*C3b*t_02_k1_4;
    w[43] = 1.0*C4b_C2a*t_03_k1_4;
    w[44] = 1.0*C4BP*PC_CRP_LF*t_04_k1_4 - 1.0*C4BP_PC_CRP_LF*t_04_k2;
}