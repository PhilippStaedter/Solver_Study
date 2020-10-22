#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Sasagawa2005(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*EGFR;
            break;
        case 1:
            dwdp[0] = 1.0*pro_EGFR;
            break;
        case 2:
            dwdp[1] = -1.0*L_EGFR;
            break;
        case 3:
            dwdp[1] = 1.0*EGF*EGFR;
            break;
        case 4:
            dwdp[2] = -1.0*L_EGFR_dimer;
            break;
        case 5:
            dwdp[2] = 1.0*pow(L_EGFR, 2);
            break;
        case 6:
            dwdp[3] = -1.0*SOS_Grb2;
            break;
        case 7:
            dwdp[3] = 1.0*Grb2*SOS;
            break;
        case 8:
            dwdp[4] = -1.0*pSOS_Grb2;
            break;
        case 9:
            dwdp[4] = 1.0*Grb2*pSOS;
            break;
        case 10:
            dwdp[5] = -1.0*L_dpEGFR;
            break;
        case 11:
            dwdp[5] = 1.0*L_EGFR_dimer;
            break;
        case 12:
            dwdp[6] = -1.0*dpEGFR_c_Cbl;
            break;
        case 13:
            dwdp[6] = 1.0*L_dpEGFR*c_Cbl;
            break;
        case 14:
            dwdp[7] = -1.0*pShc_dpEGFR;
            break;
        case 15:
            dwdp[7] = 1.0*L_dpEGFR*pShc;
            break;
        case 16:
            dwdp[8] = -1.0*Dok;
            break;
        case 17:
            dwdp[8] = 1.0*pDok;
            break;
        case 18:
            dwdp[9] = -1.0*pShc_dpEGFR_c_Cbl;
            break;
        case 19:
            dwdp[9] = 1.0*c_Cbl*pShc_dpEGFR;
            break;
        case 20:
            dwdp[10] = 1.0*pSOS;
            break;
        case 21:
            dwdp[11] = 1.0*pSOS_Grb2;
            break;
        case 22:
            dwdp[12] = -1.0*Shc_dpEGFR;
            break;
        case 23:
            dwdp[12] = 1.0*L_dpEGFR*Shc;
            break;
        case 24:
            dwdp[13] = 1.0*Shc_dpEGFR;
            break;
        case 25:
            dwdp[14] = 1.0*dpEGFR_c_Cbl;
            break;
        case 26:
            dwdp[15] = 1.0*dpEGFR_c_Cbl_ubiq;
            break;
        case 27:
            dwdp[16] = -1.0*Shc_dpEGFR_c_Cbl;
            break;
        case 28:
            dwdp[16] = 1.0*Shc_dpEGFR*c_Cbl;
            break;
        case 29:
            dwdp[17] = 1.0*Shc_dpEGFR_c_Cbl;
            break;
        case 30:
            dwdp[18] = 1.0*Shc_dpEGFR_c_Cbl_ubiq;
            break;
        case 31:
            dwdp[19] = 1.0*pShc_dpEGFR_c_Cbl;
            break;
        case 32:
            dwdp[20] = 1.0*pShc_dpEGFR_c_Cbl_ubiq;
            break;
        case 33:
            dwdp[21] = 1.0*Shc_dpEGFR_c_Cbl;
            break;
        case 34:
            dwdp[22] = -1.0*Grb2_SOS_pShc;
            break;
        case 35:
            dwdp[22] = 1.0*SOS_Grb2*pShc;
            break;
        case 36:
            dwdp[23] = -1.0*Grb2_SOS_pShc_dpEGFR;
            break;
        case 37:
            dwdp[23] = 1.0*Grb2_SOS_pShc*L_dpEGFR;
            break;
        case 38:
            dwdp[24] = -1.0*Grb2_SOS_pShc_dpEGFR;
            break;
        case 39:
            dwdp[24] = 1.0*SOS_Grb2*pShc_dpEGFR;
            break;
        case 40:
            dwdp[25] = -1.0*Grb2_SOS_pShc_dpEGFR_c_Cbl;
            break;
        case 41:
            dwdp[25] = 1.0*Grb2_SOS_pShc_dpEGFR*c_Cbl;
            break;
        case 42:
            dwdp[26] = -1.0*Grb2_SOS_pShc_dpEGFR_c_Cbl;
            break;
        case 43:
            dwdp[26] = 1.0*Grb2_SOS_pShc*dpEGFR_c_Cbl;
            break;
        case 44:
            dwdp[27] = 1.0*Grb2_SOS_pShc_dpEGFR_c_Cbl;
            break;
        case 45:
            dwdp[28] = 1.0*Grb2_SOS_pShc_dpEGFR_c_Cbl_ubiq;
            break;
        case 46:
            dwdp[29] = Grb2_SOS_pShc;
            break;
        case 47:
            dwdp[30] = -1.0*Dok*J31_J31_Vmax*(Crk_C3G_pFRS2_dpEGFR + Crk_C3G_pFRS2_dpEGFR_c_Cbl + FRS2_dpEGFR + FRS2_dpEGFR_c_Cbl + Grb2_SOS_pShc_dpEGFR + Grb2_SOS_pShc_dpEGFR_c_Cbl + L_dpEGFR + Shc_dpEGFR + Shc_dpEGFR_c_Cbl + dpEGFR_c_Cbl + pFRS2_dpEGFR + pShc_dpEGFR + pShc_dpEGFR_c_Cbl)/pow(Dok + J31_J31_Km1, 2);
            break;
        case 48:
            dwdp[30] = 1.0*Dok*(Crk_C3G_pFRS2_dpEGFR + Crk_C3G_pFRS2_dpEGFR_c_Cbl + FRS2_dpEGFR + FRS2_dpEGFR_c_Cbl + Grb2_SOS_pShc_dpEGFR + Grb2_SOS_pShc_dpEGFR_c_Cbl + L_dpEGFR + Shc_dpEGFR + Shc_dpEGFR_c_Cbl + dpEGFR_c_Cbl + pFRS2_dpEGFR + pShc_dpEGFR + pShc_dpEGFR_c_Cbl)/(Dok + J31_J31_Km1);
            break;
        case 49:
            dwdp[31] = 1.0*pShc;
            break;
        case 50:
            dwdp[32] = 1.0*pFRS2;
            break;
        case 51:
            dwdp[33] = -1.0*Crk_C3G;
            break;
        case 52:
            dwdp[33] = 1.0*C3G*Crk;
            break;
        case 53:
            dwdp[34] = -1.0*FRS2_dpEGFR;
            break;
        case 54:
            dwdp[34] = 1.0*FRS2*L_dpEGFR;
            break;
        case 55:
            dwdp[35] = -1.0*pFRS2_dpEGFR;
            break;
        case 56:
            dwdp[35] = 1.0*L_dpEGFR*pFRS2;
            break;
        case 57:
            dwdp[36] = 1.0*FRS2_dpEGFR;
            break;
        case 58:
            dwdp[37] = -1.0*Crk_C3G_pFRS2_dpEGFR;
            break;
        case 59:
            dwdp[37] = 1.0*Crk_C3G*pFRS2_dpEGFR;
            break;
        case 60:
            dwdp[38] = -1.0*FRS2_dpEGFR_c_Cbl;
            break;
        case 61:
            dwdp[38] = 1.0*FRS2_dpEGFR*c_Cbl;
            break;
        case 62:
            dwdp[39] = -1.0*pFRS2_dpEGFR_c_Cbl;
            break;
        case 63:
            dwdp[39] = 1.0*c_Cbl*pFRS2_dpEGFR;
            break;
        case 64:
            dwdp[40] = 1.0*pFRS2_dpEGFR_c_Cbl;
            break;
        case 65:
            dwdp[41] = 1.0*FRS2_dpEGFR_c_Cbl;
            break;
        case 66:
            dwdp[42] = 1.0*FRS2_dpEGFR_c_Cbl;
            break;
        case 67:
            dwdp[43] = -1.0*Crk_C3G_pFRS2_dpEGFR_c_Cbl;
            break;
        case 68:
            dwdp[43] = 1.0*Crk_C3G*pFRS2_dpEGFR_c_Cbl;
            break;
        case 69:
            dwdp[44] = 1.0*Crk_C3G_pFRS2_dpEGFR_c_Cbl;
            break;
        case 70:
            dwdp[45] = 1.0*FRS2_dpEGFR_c_Cbl_ubiq;
            break;
        case 71:
            dwdp[46] = 1.0*pFRS2_dpEGFR_c_Cbl_ubiq;
            break;
        case 72:
            dwdp[47] = -1.0*pDok_RasGAP;
            break;
        case 73:
            dwdp[47] = 1.0*RasGAP*pDok;
            break;
        case 74:
            dwdp[48] = -1.0*J50_J50_Vmax*SOS_Grb2*dppERK/pow(J50_J50_Km1 + SOS_Grb2, 2);
            break;
        case 75:
            dwdp[48] = 1.0*SOS_Grb2*dppERK/(J50_J50_Km1 + SOS_Grb2);
            break;
        case 76:
            dwdp[49] = -1.0*J51_J51_Vmax*SOS*dppERK/pow(J51_J51_Km1 + SOS, 2);
            break;
        case 77:
            dwdp[49] = 1.0*SOS*dppERK/(J51_J51_Km1 + SOS);
            break;
        case 78:
            dwdp[50] = -1.0*c_Raf_Ras_GTP;
            break;
        case 79:
            dwdp[50] = 1.0*Ras_GTP*c_Raf;
            break;
        case 80:
            dwdp[51] = -1.0*B_Raf_Rap1_GTP;
            break;
        case 81:
            dwdp[51] = 1.0*B_Raf*Rap1_GTP;
            break;
        case 82:
            dwdp[52] = -1.0*B_Raf_Ras_GTP;
            break;
        case 83:
            dwdp[52] = 1.0*B_Raf*Ras_GTP;
            break;
        case 84:
            dwdp[53] = -J57_J57_Vmax*PP2A*ppMEK/pow(J57_J57_Km1 + ppMEK, 2);
            break;
        case 85:
            dwdp[53] = PP2A*ppMEK/(J57_J57_Km1 + ppMEK);
            break;
        case 86:
            dwdp[54] = -1.0*J58_J58_Vmax*PP2A*pMEK/pow(J58_J58_Km1 + pMEK, 2);
            break;
        case 87:
            dwdp[54] = 1.0*PP2A*pMEK/(J58_J58_Km1 + pMEK);
            break;
        case 88:
            dwdp[55] = -1.0*J61_J61_Vmax*PP2A*ppMEK_ERK/pow(J61_J61_Km1 + ppMEK_ERK, 2);
            break;
        case 89:
            dwdp[55] = 1.0*PP2A*ppMEK_ERK/(J61_J61_Km1 + ppMEK_ERK);
            break;
        case 90:
            dwdp[56] = -1.0*J62_J62_Vmax*PP2A*pMEK_ERK/pow(J62_J62_Km1 + pMEK_ERK, 2);
            break;
        case 91:
            dwdp[56] = 1.0*PP2A*pMEK_ERK/(J62_J62_Km1 + pMEK_ERK);
            break;
        case 92:
            dwdp[57] = -1.0*dppERK;
            break;
        case 93:
            dwdp[57] = 1.0*pow(ppERK, 2);
            break;
        case 94:
            dwdp[58] = 1.0*Ras_GTP;
            break;
        case 95:
            dwdp[59] = 1.0*Rap1_GTP;
            break;
        case 96:
            dwdp[60] = -1.0*J68_J68_Vmax*Rap1_GDP*(Crk_C3G_pFRS2_dpEGFR + Crk_C3G_pFRS2_dpEGFR_c_Cbl + Crk_C3G_pFRS2_pTrkA_endo)/pow(J68_J68_Km1 + Rap1_GDP, 2);
            break;
        case 97:
            dwdp[60] = 1.0*Rap1_GDP*(Crk_C3G_pFRS2_dpEGFR + Crk_C3G_pFRS2_dpEGFR_c_Cbl + Crk_C3G_pFRS2_pTrkA_endo)/(J68_J68_Km1 + Rap1_GDP);
            break;
        case 98:
            dwdp[61] = -1.0*J69_J69_Vmax*Ras_GDP*(Grb2_SOS_pShc_dpEGFR + Grb2_SOS_pShc_dpEGFR_c_Cbl + Grb2_SOS_pShc_pTrkA)/pow(J69_J69_Km1 + Ras_GDP, 2);
            break;
        case 99:
            dwdp[61] = 1.0*Ras_GDP*(Grb2_SOS_pShc_dpEGFR + Grb2_SOS_pShc_dpEGFR_c_Cbl + Grb2_SOS_pShc_pTrkA)/(J69_J69_Km1 + Ras_GDP);
            break;
        case 100:
            dwdp[62] = -1.0*L_NGFR;
            break;
        case 101:
            dwdp[62] = 1.0*NGF*NGFR;
            break;
        case 102:
            dwdp[63] = 1.0*L_NGFR;
            break;
        case 103:
            dwdp[64] = 1.0*pTrkA;
            break;
        case 104:
            dwdp[65] = 1.0*pTrkA_endo;
            break;
        case 105:
            dwdp[66] = 1.0*pTrkA;
            break;
        case 106:
            dwdp[67] = -1.0*Shc_pTrkA;
            break;
        case 107:
            dwdp[67] = 1.0*Shc*pTrkA;
            break;
        case 108:
            dwdp[68] = -1.0*pShc_pTrkA;
            break;
        case 109:
            dwdp[68] = 1.0*pShc*pTrkA;
            break;
        case 110:
            dwdp[69] = -1.0*FRS2_pTrkA;
            break;
        case 111:
            dwdp[69] = 1.0*FRS2*pTrkA;
            break;
        case 112:
            dwdp[70] = -1.0*pFRS2_pTrkA;
            break;
        case 113:
            dwdp[70] = 1.0*pFRS2*pTrkA;
            break;
        case 114:
            dwdp[71] = -1.0*Shc_pTrkA_endo;
            break;
        case 115:
            dwdp[71] = 1.0*Shc*pTrkA_endo;
            break;
        case 116:
            dwdp[72] = -1.0*pShc_pTrkA_endo;
            break;
        case 117:
            dwdp[72] = 1.0*pShc*pTrkA_endo;
            break;
        case 118:
            dwdp[73] = 1.0*Shc_pTrkA_endo;
            break;
        case 119:
            dwdp[74] = 1.0*Shc_pTrkA;
            break;
        case 120:
            dwdp[75] = 1.0*FRS2_pTrkA;
            break;
        case 121:
            dwdp[76] = -1.0*FRS2_pTrkA_endo;
            break;
        case 122:
            dwdp[76] = 1.0*FRS2*pTrkA_endo;
            break;
        case 123:
            dwdp[77] = -1.0*pFRS2_pTrkA_endo;
            break;
        case 124:
            dwdp[77] = 1.0*pFRS2*pTrkA_endo;
            break;
        case 125:
            dwdp[78] = 1.0*FRS2_pTrkA_endo;
            break;
        case 126:
            dwdp[79] = 1.0*FRS2_pTrkA;
            break;
        case 127:
            dwdp[80] = 1.0*pFRS2_pTrkA;
            break;
        case 128:
            dwdp[81] = 1.0*Shc_pTrkA;
            break;
        case 129:
            dwdp[82] = 1.0*pShc_pTrkA;
            break;
        case 130:
            dwdp[83] = 1.0*FRS2_pTrkA_endo;
            break;
        case 131:
            dwdp[84] = 1.0*Shc_pTrkA_endo;
            break;
        case 132:
            dwdp[85] = 1.0*pShc_pTrkA_endo;
            break;
        case 133:
            dwdp[86] = -1.0*Grb2_SOS_pShc_pTrkA_endo;
            break;
        case 134:
            dwdp[86] = 1.0*SOS_Grb2*pShc_pTrkA_endo;
            break;
        case 135:
            dwdp[87] = -1.0*Grb2_SOS_pShc_pTrkA;
            break;
        case 136:
            dwdp[87] = 1.0*SOS_Grb2*pShc_pTrkA;
            break;
        case 137:
            dwdp[88] = 1.0*Grb2_SOS_pShc_pTrkA;
            break;
        case 138:
            dwdp[89] = 1.0*Crk_C3G_pFRS2_pTrkA;
            break;
        case 139:
            dwdp[90] = 1.0*pFRS2_pTrkA;
            break;
        case 140:
            dwdp[91] = 1.0*FRS2_pTrkA;
            break;
        case 141:
            dwdp[92] = 1.0*pShc_pTrkA;
            break;
        case 142:
            dwdp[93] = 1.0*Shc_pTrkA;
            break;
        case 143:
            dwdp[94] = -1.0*Crk_C3G_pFRS2_pTrkA;
            break;
        case 144:
            dwdp[94] = 1.0*Crk_C3G*pFRS2_pTrkA;
            break;
        case 145:
            dwdp[95] = -1.0*Crk_C3G_pFRS2_pTrkA_endo;
            break;
        case 146:
            dwdp[95] = 1.0*Crk_C3G*pFRS2_pTrkA_endo;
            break;
        case 147:
            dwdp[96] = -1.0*Grb2_SOS_pShc_pTrkA;
            break;
        case 148:
            dwdp[96] = 1.0*Grb2_SOS_pShc*pTrkA;
            break;
        case 149:
            dwdp[97] = -1.0*Grb2_SOS_pShc_pTrkA_endo;
            break;
        case 150:
            dwdp[97] = 1.0*Grb2_SOS_pShc*pTrkA_endo;
            break;
        case 151:
            dwdp[98] = 1.0*Crk_C3G_pFRS2_pTrkA;
            break;
        case 152:
            dwdp[99] = 1.0*Crk_C3G_pFRS2_pTrkA_endo;
            break;
        case 153:
            dwdp[100] = 1.0*Grb2_SOS_pShc_pTrkA;
            break;
        case 154:
            dwdp[101] = 1.0*Grb2_SOS_pShc_pTrkA_endo;
            break;
        case 155:
            dwdp[102] = pFRS2_pTrkA_endo;
            break;
        case 156:
            dwdp[103] = -1.0*NGFR;
            break;
        case 157:
            dwdp[103] = 1.0*pro_TrkA;
            break;
        case 158:
            dwdp[104] = -1.0*Shc_dpEGFR_c_Cbl;
            break;
        case 159:
            dwdp[104] = 1.0*Shc*dpEGFR_c_Cbl;
            break;
        case 160:
            dwdp[105] = -1.0*pShc_dpEGFR_c_Cbl;
            break;
        case 161:
            dwdp[105] = 1.0*dpEGFR_c_Cbl*pShc;
            break;
        case 162:
            dwdp[106] = -1.0*Grb2_SOS_pShc_dpEGFR_c_Cbl;
            break;
        case 163:
            dwdp[106] = 1.0*SOS_Grb2*pShc_dpEGFR_c_Cbl;
            break;
        case 164:
            dwdp[107] = -1.0*Crk_C3G_pFRS2_dpEGFR_c_Cbl;
            break;
        case 165:
            dwdp[107] = 1.0*Crk_C3G_pFRS2_dpEGFR*c_Cbl;
            break;
        case 166:
            dwdp[108] = -1.0*FRS2_dpEGFR_c_Cbl;
            break;
        case 167:
            dwdp[108] = 1.0*FRS2*dpEGFR_c_Cbl;
            break;
        case 168:
            dwdp[109] = -1.0*pFRS2_dpEGFR_c_Cbl;
            break;
        case 169:
            dwdp[109] = 1.0*dpEGFR_c_Cbl*pFRS2;
            break;
        case 170:
            dwdp[110] = -1.0*J121_J121_Vmax*Ras_GTP*pDok_RasGAP/pow(J121_J121_Km1 + Ras_GTP, 2);
            break;
        case 171:
            dwdp[110] = 1.0*Ras_GTP*pDok_RasGAP/(J121_J121_Km1 + Ras_GTP);
            break;
        case 172:
            dwdp[111] = -1.0*J122_J122_Vmax*Rap1GAP*Rap1_GTP/pow(J122_J122_Km1 + Rap1_GTP, 2);
            break;
        case 173:
            dwdp[111] = 1.0*Rap1GAP*Rap1_GTP/(J122_J122_Km1 + Rap1_GTP);
            break;
        case 174:
            dwdp[112] = -1.0*Dok*J123_J123_Vmax*(Crk_C3G_pFRS2_pTrkA + FRS2_pTrkA + Grb2_SOS_pShc_pTrkA + Shc_pTrkA + pFRS2_pTrkA + pShc_pTrkA + pTrkA)/pow(Dok + J123_J123_Km1, 2);
            break;
        case 175:
            dwdp[112] = 1.0*Dok*(Crk_C3G_pFRS2_pTrkA + FRS2_pTrkA + Grb2_SOS_pShc_pTrkA + Shc_pTrkA + pFRS2_pTrkA + pShc_pTrkA + pTrkA)/(Dok + J123_J123_Km1);
            break;
        case 176:
            dwdp[113] = -1.0*Grb2_SOS_pShc*J124_J124_Vmax*dppERK/pow(Grb2_SOS_pShc + J124_J124_Km1, 2);
            break;
        case 177:
            dwdp[113] = 1.0*Grb2_SOS_pShc*dppERK/(Grb2_SOS_pShc + J124_J124_Km1);
            break;
        case 178:
            dwdp[114] = -1.0*MEK_ERK;
            break;
        case 179:
            dwdp[114] = 1.0*ERK*MEK;
            break;
        case 180:
            dwdp[115] = -1.0*pMEK_ERK;
            break;
        case 181:
            dwdp[115] = 1.0*ERK*pMEK;
            break;
        case 182:
            dwdp[116] = -1.0*ppMEK_ERK;
            break;
        case 183:
            dwdp[116] = 1.0*ERK*ppMEK;
            break;
        case 184:
            dwdp[117] = 1.0*ppMEK_ERK;
            break;
        case 185:
            dwdp[118] = -1.0*J137_J137_Vmax*c_Raf_Ras_GTP*pDok_RasGAP/pow(J137_J137_Km1 + c_Raf_Ras_GTP, 2);
            break;
        case 186:
            dwdp[118] = 1.0*c_Raf_Ras_GTP*pDok_RasGAP/(J137_J137_Km1 + c_Raf_Ras_GTP);
            break;
        case 187:
            dwdp[119] = -1.0*B_Raf_Ras_GTP*J138_J138_Vmax*pDok_RasGAP/pow(B_Raf_Ras_GTP + J138_J138_Km1, 2);
            break;
        case 188:
            dwdp[119] = 1.0*B_Raf_Ras_GTP*pDok_RasGAP/(B_Raf_Ras_GTP + J138_J138_Km1);
            break;
        case 189:
            dwdp[120] = -1.0*B_Raf_Rap1_GTP*J139_J139_Vmax*Rap1GAP/pow(B_Raf_Rap1_GTP + J139_J139_Km1, 2);
            break;
        case 190:
            dwdp[120] = 1.0*B_Raf_Rap1_GTP*Rap1GAP/(B_Raf_Rap1_GTP + J139_J139_Km1);
            break;
        case 191:
            dwdp[121] = -1.0*c_Raf_Ras_GTP_MEK;
            break;
        case 192:
            dwdp[121] = 1.0*MEK*c_Raf_Ras_GTP;
            break;
        case 193:
            dwdp[122] = -1.0*c_Raf_Ras_GTP_pMEK;
            break;
        case 194:
            dwdp[122] = 1.0*c_Raf_Ras_GTP*pMEK;
            break;
        case 195:
            dwdp[123] = -1.0*c_Raf_Ras_GTP_MEK_ERK;
            break;
        case 196:
            dwdp[123] = 1.0*MEK_ERK*c_Raf_Ras_GTP;
            break;
        case 197:
            dwdp[124] = -1.0*c_Raf_Ras_GTP_pMEK_ERK;
            break;
        case 198:
            dwdp[124] = 1.0*c_Raf_Ras_GTP*pMEK_ERK;
            break;
        case 199:
            dwdp[125] = -1.0*B_Raf_Ras_GTP_MEK;
            break;
        case 200:
            dwdp[125] = 1.0*B_Raf_Ras_GTP*MEK;
            break;
        case 201:
            dwdp[126] = -1.0*B_Raf_Ras_GTP_pMEK;
            break;
        case 202:
            dwdp[126] = 1.0*B_Raf_Ras_GTP*pMEK;
            break;
        case 203:
            dwdp[127] = -1.0*B_Raf_Ras_GTP_MEK_ERK;
            break;
        case 204:
            dwdp[127] = 1.0*B_Raf_Ras_GTP*MEK_ERK;
            break;
        case 205:
            dwdp[128] = -1.0*B_Raf_Ras_GTP_pMEK_ERK;
            break;
        case 206:
            dwdp[128] = 1.0*B_Raf_Ras_GTP*pMEK_ERK;
            break;
        case 207:
            dwdp[129] = -1.0*B_Raf_Rap1_GTP_MEK;
            break;
        case 208:
            dwdp[129] = 1.0*B_Raf_Rap1_GTP*MEK;
            break;
        case 209:
            dwdp[130] = -1.0*B_Raf_Rap1_GTP_pMEK;
            break;
        case 210:
            dwdp[130] = 1.0*B_Raf_Rap1_GTP*pMEK;
            break;
        case 211:
            dwdp[131] = -1.0*B_Raf_Rap1_GTP_MEK_ERK;
            break;
        case 212:
            dwdp[131] = 1.0*B_Raf_Rap1_GTP*MEK_ERK;
            break;
        case 213:
            dwdp[132] = -1.0*B_Raf_Rap1_GTP_pMEK_ERK;
            break;
        case 214:
            dwdp[132] = 1.0*B_Raf_Rap1_GTP*pMEK_ERK;
            break;
        case 215:
            dwdp[133] = 1.0*c_Raf_Ras_GTP_MEK;
            break;
        case 216:
            dwdp[134] = 1.0*c_Raf_Ras_GTP_pMEK;
            break;
        case 217:
            dwdp[135] = 1.0*c_Raf_Ras_GTP_MEK_ERK;
            break;
        case 218:
            dwdp[136] = 1.0*c_Raf_Ras_GTP_pMEK_ERK;
            break;
        case 219:
            dwdp[137] = 1.0*B_Raf_Ras_GTP_MEK;
            break;
        case 220:
            dwdp[138] = 1.0*B_Raf_Ras_GTP_pMEK;
            break;
        case 221:
            dwdp[139] = 1.0*B_Raf_Ras_GTP_MEK_ERK;
            break;
        case 222:
            dwdp[140] = 1.0*B_Raf_Ras_GTP_pMEK_ERK;
            break;
        case 223:
            dwdp[141] = 1.0*B_Raf_Rap1_GTP_MEK;
            break;
        case 224:
            dwdp[142] = 1.0*B_Raf_Rap1_GTP_pMEK;
            break;
        case 225:
            dwdp[143] = 1.0*B_Raf_Rap1_GTP_MEK_ERK;
            break;
        case 226:
            dwdp[144] = 1.0*B_Raf_Rap1_GTP_pMEK_ERK;
            break;
        case 227:
            dwdp[145] = 1.0*Crk_C3G_pFRS2_dpEGFR_c_Cbl_ubiq;
            break;
        case 228:
            dwdp[146] = -1.0*dppERK_MKP3;
            break;
        case 229:
            dwdp[146] = 1.0*MKP3*dppERK;
            break;
        case 230:
            dwdp[147] = -1.0*ppERK_MKP3;
            break;
        case 231:
            dwdp[147] = 1.0*MKP3*ppERK;
            break;
        case 232:
            dwdp[148] = 1.0*ppERK_MKP3;
            break;
        case 233:
            dwdp[149] = 1.0*dppERK_MKP3;
            break;
    }
}