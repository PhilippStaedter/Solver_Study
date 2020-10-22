#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_liu2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*PC_CRP*ka02_1;
    dwdx[1] = 1.0*PC_CRP_LF*kd05_1;
    dwdx[2] = 1.0*PC_CRP_LF_MASP*kd09_1;
    dwdx[3] = 1.0*GlcNac_LF_CRP*ke02_1;
    dwdx[4] = -1.0*C2*PC_CRP_C1*ka04_1/pow(C2 + ka04_2, 2) + 1.0*PC_CRP_C1*ka04_1/(C2 + ka04_2);
    dwdx[5] = -1.0*C2*GlcNac_LF_MASP*kb04_1/pow(C2 + kb04_2, 2) + 1.0*GlcNac_LF_MASP*kb04_1/(C2 + kb04_2);
    dwdx[6] = -1.0*C2*PC_CRP_LF_MASP*kd04_1/pow(C2 + kd04_2, 2) + 1.0*PC_CRP_LF_MASP*kd04_1/(C2 + kd04_2);
    dwdx[7] = -1.0*C2*PC_CRP_LF_C1*kd07_1/pow(C2 + kd07_2, 2) + 1.0*PC_CRP_LF_C1*kd07_1/(C2 + kd07_2);
    dwdx[8] = -1.0*C2*PC_CRP_LF_C1_MASP*kd11_1/pow(C2 + kd11_2, 2) + 1.0*PC_CRP_LF_C1_MASP*kd11_1/(C2 + kd11_2);
    dwdx[9] = -1.0*C2*GlcNac_LF_CRP_C1*ke04_1/pow(C2 + ke04_2, 2) + 1.0*GlcNac_LF_CRP_C1*ke04_1/(C2 + ke04_2);
    dwdx[10] = -1.0*C2*GlcNac_LF_CRP_MASP*ke07_1/pow(C2 + ke07_2, 2) + 1.0*GlcNac_LF_CRP_MASP*ke07_1/(C2 + ke07_2);
    dwdx[11] = -1.0*C2*GlcNac_HF_MASP*kg04_1/pow(C2 + kg04_2, 2) + 1.0*GlcNac_HF_MASP*kg04_1/(C2 + kg04_2);
    dwdx[12] = 1.0*C4b*kc01_1;
    dwdx[13] = 1.0*C4b_C2a*kc02;
    dwdx[14] = 1.0*kc03_1;
    dwdx[15] = 1.0*t_02_k1_4;
    dwdx[16] = -1.0*C4*PC_CRP_C1*ka03_1/pow(C4 + ka03_2, 2) + 1.0*PC_CRP_C1*ka03_1/(C4 + ka03_2);
    dwdx[17] = -1.0*C4*GlcNac_LF_MASP*kb03_1/pow(C4 + kb03_2, 2) + 1.0*GlcNac_LF_MASP*kb03_1/(C4 + kb03_2);
    dwdx[18] = -1.0*C4*PC_CRP_LF_MASP*kd03_1/pow(C4 + kd03_2, 2) + 1.0*PC_CRP_LF_MASP*kd03_1/(C4 + kd03_2);
    dwdx[19] = -1.0*C4*PC_CRP_LF_C1*kd06_1/pow(C4 + kd06_2, 2) + 1.0*PC_CRP_LF_C1*kd06_1/(C4 + kd06_2);
    dwdx[20] = -1.0*C4*PC_CRP_LF_C1_MASP*kd10_1/pow(C4 + kd10_2, 2) + 1.0*PC_CRP_LF_C1_MASP*kd10_1/(C4 + kd10_2);
    dwdx[21] = -1.0*C4*GlcNac_LF_CRP_C1*ke03_1/pow(C4 + ke03_2, 2) + 1.0*GlcNac_LF_CRP_C1*ke03_1/(C4 + ke03_2);
    dwdx[22] = -1.0*C4*GlcNac_LF_CRP_MASP*ke06_1/pow(C4 + ke06_2, 2) + 1.0*GlcNac_LF_CRP_MASP*ke06_1/(C4 + ke06_2);
    dwdx[23] = -1.0*C4*GlcNac_HF_MASP*kg03_1/pow(C4 + kg03_2, 2) + 1.0*GlcNac_HF_MASP*kg03_1/(C4 + kg03_2);
    dwdx[24] = 1.0*PC_CRP*kf01_1;
    dwdx[25] = 1.0*GlcNac_LF_CRP*kf02_1;
    dwdx[26] = 1.0*C4b_C2a*kf03;
    dwdx[27] = 1.0*C4b*kf04_1;
    dwdx[28] = 1.0*C4b_C2a*kf05;
    dwdx[29] = 1.0*C4b_C2a*kf06_1;
    dwdx[30] = 1.0*dC4b_C2a*kf07_1;
    dwdx[31] = 1.0*t_01_k1_4;
    dwdx[32] = 1.0*PC_CRP_LF*t_04_k1_4;
    dwdx[33] = -1.0*kf04_2;
    dwdx[34] = -1.0*kf02_2;
    dwdx[35] = -1.0*kf01_2;
    dwdx[36] = -1.0*t_04_k2;
    dwdx[37] = 1.0*C2a*kc01_1;
    dwdx[38] = 1.0*C4BP*kf04_1;
    dwdx[39] = -1.0*kc01_2;
    dwdx[40] = 1.0*C3*kc02;
    dwdx[41] = 1.0*kc04_1;
    dwdx[42] = 1.0*C4BP*kf03;
    dwdx[43] = 1.0*C4BP*kf05;
    dwdx[44] = 1.0*C4BP*kf06_1;
    dwdx[45] = 1.0*t_03_k1_4;
    dwdx[46] = -1.0*kf06_2;
    dwdx[47] = 1.0*PC*ka01_1;
    dwdx[48] = 1.0*GlcNac_LF*ke01_1;
    dwdx[49] = 1.0*LF*kb01_1;
    dwdx[50] = -1.0*kg01_2;
    dwdx[51] = 1.0*MASP*kg02_1;
    dwdx[52] = -1.0*kg02_2;
    dwdx[53] = 1.0*C4*kg03_1/(C4 + kg03_2);
    dwdx[54] = 1.0*C2*kg04_1/(C2 + kg04_2);
    dwdx[55] = -1.0*kb01_2;
    dwdx[56] = 1.0*MASP*kb02_1;
    dwdx[57] = 1.0*CRP*ke01_1;
    dwdx[58] = -1.0*ke01_2;
    dwdx[59] = 1.0*C1*ke02_1;
    dwdx[60] = 1.0*MASP*ke05_1;
    dwdx[61] = 1.0*C4BP*kf02_1;
    dwdx[62] = -1.0*ke02_2;
    dwdx[63] = 1.0*C4*ke03_1/(C4 + ke03_2);
    dwdx[64] = 1.0*C2*ke04_1/(C2 + ke04_2);
    dwdx[65] = -1.0*ke05_2;
    dwdx[66] = 1.0*C4*ke06_1/(C4 + ke06_2);
    dwdx[67] = 1.0*C2*ke07_1/(C2 + ke07_2);
    dwdx[68] = -1.0*kb02_2;
    dwdx[69] = 1.0*C4*kb03_1/(C4 + kb03_2);
    dwdx[70] = 1.0*C2*kb04_1/(C2 + kb04_2);
    dwdx[71] = 1.0*X*kg01_1;
    dwdx[72] = 1.0*GlcNac*kb01_1;
    dwdx[73] = 1.0*PC_CRP*kd01_1;
    dwdx[74] = 1.0*GlcNac_LF*kb02_1;
    dwdx[75] = 1.0*PC_CRP_LF*kd02_1;
    dwdx[76] = 1.0*PC_CRP_LF_C1*kd08_1;
    dwdx[77] = 1.0*GlcNac_LF_CRP*ke05_1;
    dwdx[78] = 1.0*GlcNac_HF*kg02_1;
    dwdx[79] = 1.0*CRP*ka01_1;
    dwdx[80] = -1.0*ka01_2;
    dwdx[81] = 1.0*C1*ka02_1;
    dwdx[82] = 1.0*LF*kd01_1;
    dwdx[83] = 1.0*C4BP*kf01_1;
    dwdx[84] = -1.0*ka02_2;
    dwdx[85] = 1.0*C4*ka03_1/(C4 + ka03_2);
    dwdx[86] = 1.0*C2*ka04_1/(C2 + ka04_2);
    dwdx[87] = -1.0*kd01_2;
    dwdx[88] = 1.0*MASP*kd02_1;
    dwdx[89] = 1.0*C1*kd05_1;
    dwdx[90] = 1.0*C4BP*t_04_k1_4;
    dwdx[91] = -1.0*kd05_2;
    dwdx[92] = 1.0*C4*kd06_1/(C4 + kd06_2);
    dwdx[93] = 1.0*C2*kd07_1/(C2 + kd07_2);
    dwdx[94] = 1.0*MASP*kd08_1;
    dwdx[95] = -1.0*kd08_2;
    dwdx[96] = -1.0*kd09_2;
    dwdx[97] = 1.0*C4*kd10_1/(C4 + kd10_2);
    dwdx[98] = 1.0*C2*kd11_1/(C2 + kd11_2);
    dwdx[99] = -1.0*kd02_2;
    dwdx[100] = 1.0*C4*kd03_1/(C4 + kd03_2);
    dwdx[101] = 1.0*C2*kd04_1/(C2 + kd04_2);
    dwdx[102] = 1.0*C1*kd09_1;
    dwdx[103] = 1.0*HF*kg01_1;
    dwdx[104] = -1.0*kc03_2;
    dwdx[105] = -1.0*kc04_2;
    dwdx[106] = 1.0*C4BP*kf07_1;
    dwdx[107] = -1.0*kf07_2;
}