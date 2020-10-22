#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_liu2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*CRP*PC;
            break;
        case 1:
            dwdp[0] = -1.0*PC_CRP;
            break;
        case 2:
            dwdp[1] = 1.0*C1*PC_CRP;
            break;
        case 3:
            dwdp[1] = -1.0*PC_CRP_C1;
            break;
        case 4:
            dwdp[2] = 1.0*C4*PC_CRP_C1/(C4 + ka03_2);
            break;
        case 5:
            dwdp[2] = -1.0*C4*PC_CRP_C1*ka03_1/pow(C4 + ka03_2, 2);
            break;
        case 6:
            dwdp[3] = 1.0*C2*PC_CRP_C1/(C2 + ka04_2);
            break;
        case 7:
            dwdp[3] = -1.0*C2*PC_CRP_C1*ka04_1/pow(C2 + ka04_2, 2);
            break;
        case 8:
            dwdp[4] = 1.0*GlcNac*LF;
            break;
        case 9:
            dwdp[4] = -1.0*GlcNac_LF;
            break;
        case 10:
            dwdp[5] = 1.0*GlcNac_LF*MASP;
            break;
        case 11:
            dwdp[5] = -1.0*GlcNac_LF_MASP;
            break;
        case 12:
            dwdp[6] = 1.0*C4*GlcNac_LF_MASP/(C4 + kb03_2);
            break;
        case 13:
            dwdp[6] = -1.0*C4*GlcNac_LF_MASP*kb03_1/pow(C4 + kb03_2, 2);
            break;
        case 14:
            dwdp[7] = 1.0*C2*GlcNac_LF_MASP/(C2 + kb04_2);
            break;
        case 15:
            dwdp[7] = -1.0*C2*GlcNac_LF_MASP*kb04_1/pow(C2 + kb04_2, 2);
            break;
        case 16:
            dwdp[8] = 1.0*C2a*C4b;
            break;
        case 17:
            dwdp[8] = -1.0*C4b_C2a;
            break;
        case 18:
            dwdp[9] = 1.0*C3*C4b_C2a;
            break;
        case 19:
            dwdp[10] = 1.0*C3b;
            break;
        case 20:
            dwdp[10] = -1.0*dC3b;
            break;
        case 21:
            dwdp[11] = 1.0*C4b_C2a;
            break;
        case 22:
            dwdp[11] = -1.0*dC4b_C2a;
            break;
        case 23:
            dwdp[12] = 1.0*LF*PC_CRP;
            break;
        case 24:
            dwdp[12] = -1.0*PC_CRP_LF;
            break;
        case 25:
            dwdp[13] = 1.0*MASP*PC_CRP_LF;
            break;
        case 26:
            dwdp[13] = -1.0*PC_CRP_LF_MASP;
            break;
        case 27:
            dwdp[14] = 1.0*C4*PC_CRP_LF_MASP/(C4 + kd03_2);
            break;
        case 28:
            dwdp[14] = -1.0*C4*PC_CRP_LF_MASP*kd03_1/pow(C4 + kd03_2, 2);
            break;
        case 29:
            dwdp[15] = 1.0*C2*PC_CRP_LF_MASP/(C2 + kd04_2);
            break;
        case 30:
            dwdp[15] = -1.0*C2*PC_CRP_LF_MASP*kd04_1/pow(C2 + kd04_2, 2);
            break;
        case 31:
            dwdp[16] = 1.0*C1*PC_CRP_LF;
            break;
        case 32:
            dwdp[16] = -1.0*PC_CRP_LF_C1;
            break;
        case 33:
            dwdp[17] = 1.0*C4*PC_CRP_LF_C1/(C4 + kd06_2);
            break;
        case 34:
            dwdp[17] = -1.0*C4*PC_CRP_LF_C1*kd06_1/pow(C4 + kd06_2, 2);
            break;
        case 35:
            dwdp[18] = 1.0*C2*PC_CRP_LF_C1/(C2 + kd07_2);
            break;
        case 36:
            dwdp[18] = -1.0*C2*PC_CRP_LF_C1*kd07_1/pow(C2 + kd07_2, 2);
            break;
        case 37:
            dwdp[19] = 1.0*MASP*PC_CRP_LF_C1;
            break;
        case 38:
            dwdp[19] = -1.0*PC_CRP_LF_C1_MASP;
            break;
        case 39:
            dwdp[20] = 1.0*C1*PC_CRP_LF_MASP;
            break;
        case 40:
            dwdp[20] = -1.0*PC_CRP_LF_C1_MASP;
            break;
        case 41:
            dwdp[21] = 1.0*C4*PC_CRP_LF_C1_MASP/(C4 + kd10_2);
            break;
        case 42:
            dwdp[21] = -1.0*C4*PC_CRP_LF_C1_MASP*kd10_1/pow(C4 + kd10_2, 2);
            break;
        case 43:
            dwdp[22] = 1.0*C2*PC_CRP_LF_C1_MASP/(C2 + kd11_2);
            break;
        case 44:
            dwdp[22] = -1.0*C2*PC_CRP_LF_C1_MASP*kd11_1/pow(C2 + kd11_2, 2);
            break;
        case 45:
            dwdp[23] = 1.0*CRP*GlcNac_LF;
            break;
        case 46:
            dwdp[23] = -1.0*GlcNac_LF_CRP;
            break;
        case 47:
            dwdp[24] = 1.0*C1*GlcNac_LF_CRP;
            break;
        case 48:
            dwdp[24] = -1.0*GlcNac_LF_CRP_C1;
            break;
        case 49:
            dwdp[25] = 1.0*C4*GlcNac_LF_CRP_C1/(C4 + ke03_2);
            break;
        case 50:
            dwdp[25] = -1.0*C4*GlcNac_LF_CRP_C1*ke03_1/pow(C4 + ke03_2, 2);
            break;
        case 51:
            dwdp[26] = 1.0*C2*GlcNac_LF_CRP_C1/(C2 + ke04_2);
            break;
        case 52:
            dwdp[26] = -1.0*C2*GlcNac_LF_CRP_C1*ke04_1/pow(C2 + ke04_2, 2);
            break;
        case 53:
            dwdp[27] = 1.0*GlcNac_LF_CRP*MASP;
            break;
        case 54:
            dwdp[27] = -1.0*GlcNac_LF_CRP_MASP;
            break;
        case 55:
            dwdp[28] = 1.0*C4*GlcNac_LF_CRP_MASP/(C4 + ke06_2);
            break;
        case 56:
            dwdp[28] = -1.0*C4*GlcNac_LF_CRP_MASP*ke06_1/pow(C4 + ke06_2, 2);
            break;
        case 57:
            dwdp[29] = 1.0*C2*GlcNac_LF_CRP_MASP/(C2 + ke07_2);
            break;
        case 58:
            dwdp[29] = -1.0*C2*GlcNac_LF_CRP_MASP*ke07_1/pow(C2 + ke07_2, 2);
            break;
        case 59:
            dwdp[30] = 1.0*C4BP*PC_CRP;
            break;
        case 60:
            dwdp[30] = -1.0*C4BP_PC_CRP;
            break;
        case 61:
            dwdp[31] = 1.0*C4BP*GlcNac_LF_CRP;
            break;
        case 62:
            dwdp[31] = -1.0*C4BP_GlcNac_LF_CRP;
            break;
        case 63:
            dwdp[32] = 1.0*C4BP*C4b_C2a;
            break;
        case 64:
            dwdp[33] = 1.0*C4BP*C4b;
            break;
        case 65:
            dwdp[33] = -1.0*C4BP_C4b;
            break;
        case 66:
            dwdp[34] = 1.0*C4BP*C4b_C2a;
            break;
        case 67:
            dwdp[35] = 1.0*C4BP*C4b_C2a;
            break;
        case 68:
            dwdp[35] = -1.0*C4b_C2a_C4BP;
            break;
        case 69:
            dwdp[36] = 1.0*C4BP*dC4b_C2a;
            break;
        case 70:
            dwdp[36] = -1.0*dC4b_C2a_C4BP;
            break;
        case 71:
            dwdp[37] = 1.0*HF*X;
            break;
        case 72:
            dwdp[37] = -1.0*GlcNac_HF;
            break;
        case 73:
            dwdp[38] = 1.0*GlcNac_HF*MASP;
            break;
        case 74:
            dwdp[38] = -1.0*GlcNac_HF_MASP;
            break;
        case 75:
            dwdp[39] = 1.0*C4*GlcNac_HF_MASP/(C4 + kg03_2);
            break;
        case 76:
            dwdp[39] = -1.0*C4*GlcNac_HF_MASP*kg03_1/pow(C4 + kg03_2, 2);
            break;
        case 77:
            dwdp[40] = 1.0*C2*GlcNac_HF_MASP/(C2 + kg04_2);
            break;
        case 78:
            dwdp[40] = -1.0*C2*GlcNac_HF_MASP*kg04_1/pow(C2 + kg04_2, 2);
            break;
        case 84:
            dwdp[41] = 1.0*C4BP;
            break;
        case 85:
            dwdp[42] = 1.0*C3b;
            break;
        case 86:
            dwdp[43] = 1.0*C4b_C2a;
            break;
        case 87:
            dwdp[44] = -1.0*C4BP_PC_CRP_LF;
            break;
        case 88:
            dwdp[44] = 1.0*C4BP*PC_CRP_LF;
            break;
    }
}