#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Sasagawa2005(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EGFR;
    x_rdata[1] = L_EGFR;
    x_rdata[2] = L_EGFR_dimer;
    x_rdata[3] = SOS;
    x_rdata[4] = L_dpEGFR;
    x_rdata[5] = pSOS;
    x_rdata[6] = SOS_Grb2;
    x_rdata[7] = Grb2;
    x_rdata[8] = Dok;
    x_rdata[9] = pDok;
    x_rdata[10] = Crk;
    x_rdata[11] = FRS2;
    x_rdata[12] = Shc;
    x_rdata[13] = pSOS_Grb2;
    x_rdata[14] = Rap1_GDP;
    x_rdata[15] = MEK;
    x_rdata[16] = MKP3;
    x_rdata[17] = pShc_dpEGFR;
    x_rdata[18] = dpEGFR_c_Cbl;
    x_rdata[19] = B_Raf_Rap1_GTP;
    x_rdata[20] = pShc_dpEGFR_c_Cbl;
    x_rdata[21] = pFRS2_dpEGFR_c_Cbl;
    x_rdata[22] = Shc_dpEGFR;
    x_rdata[23] = c_Cbl;
    x_rdata[24] = RasGAP;
    x_rdata[25] = c_Raf;
    x_rdata[26] = B_Raf;
    x_rdata[27] = ERK;
    x_rdata[28] = PP2A;
    x_rdata[29] = Ras_GDP;
    x_rdata[30] = Rap1GAP;
    x_rdata[31] = C3G;
    x_rdata[32] = NGFR;
    x_rdata[33] = pShc;
    x_rdata[34] = pFRS2_dpEGFR;
    x_rdata[35] = pTrkA_endo;
    x_rdata[36] = MEK_ERK;
    x_rdata[37] = pMEK_ERK;
    x_rdata[38] = FRS2_dpEGFR_c_Cbl_ubiq;
    x_rdata[39] = Crk_C3G_pFRS2_dpEGFR_c_Cbl;
    x_rdata[40] = pShc_dpEGFR_c_Cbl_ubiq;
    x_rdata[41] = Crk_C3G_pFRS2_dpEGFR;
    x_rdata[42] = Grb2_SOS_pShc_dpEGFR_c_Cbl_ubiq;
    x_rdata[43] = Grb2_SOS_pShc_dpEGFR_c_Cbl;
    x_rdata[44] = Shc_dpEGFR_c_Cbl_ubiq;
    x_rdata[45] = dpEGFR_c_Cbl_ubiq;
    x_rdata[46] = proteosome;
    x_rdata[47] = Grb2_SOS_pShc;
    x_rdata[48] = Shc_dpEGFR_c_Cbl;
    x_rdata[49] = Grb2_SOS_pShc_dpEGFR;
    x_rdata[50] = pFRS2;
    x_rdata[51] = FRS2_dpEGFR;
    x_rdata[52] = pDok_RasGAP;
    x_rdata[53] = pMEK;
    x_rdata[54] = FRS2_dpEGFR_c_Cbl;
    x_rdata[55] = pFRS2_dpEGFR_c_Cbl_ubiq;
    x_rdata[56] = Ras_GTP;
    x_rdata[57] = Crk_C3G_pFRS2_dpEGFR_c_Cbl_ubiq;
    x_rdata[58] = c_Raf_Ras_GTP;
    x_rdata[59] = B_Raf_Ras_GTP;
    x_rdata[60] = ppMEK;
    x_rdata[61] = ppERK;
    x_rdata[62] = pTrkA;
    x_rdata[63] = Crk_C3G;
    x_rdata[64] = Rap1_GTP;
    x_rdata[65] = L_NGFR;
    x_rdata[66] = ppMEK_ERK;
    x_rdata[67] = dppERK;
    x_rdata[68] = Shc_pTrkA;
    x_rdata[69] = Shc_pTrkA_endo;
    x_rdata[70] = pShc_pTrkA;
    x_rdata[71] = pFRS2_pTrkA;
    x_rdata[72] = FRS2_pTrkA;
    x_rdata[73] = pShc_pTrkA_endo;
    x_rdata[74] = FRS2_pTrkA_endo;
    x_rdata[75] = pFRS2_pTrkA_endo;
    x_rdata[76] = Crk_C3G_pFRS2_pTrkA_endo;
    x_rdata[77] = Grb2_SOS_pShc_pTrkA;
    x_rdata[78] = Crk_C3G_pFRS2_pTrkA;
    x_rdata[79] = Grb2_SOS_pShc_pTrkA_endo;
    x_rdata[80] = c_Raf_Ras_GTP_MEK;
    x_rdata[81] = c_Raf_Ras_GTP_pMEK;
    x_rdata[82] = c_Raf_Ras_GTP_MEK_ERK;
    x_rdata[83] = c_Raf_Ras_GTP_pMEK_ERK;
    x_rdata[84] = B_Raf_Ras_GTP_MEK;
    x_rdata[85] = B_Raf_Ras_GTP_pMEK;
    x_rdata[86] = B_Raf_Ras_GTP_MEK_ERK;
    x_rdata[87] = B_Raf_Ras_GTP_pMEK_ERK;
    x_rdata[88] = B_Raf_Rap1_GTP_MEK;
    x_rdata[89] = B_Raf_Rap1_GTP_pMEK;
    x_rdata[90] = B_Raf_Rap1_GTP_MEK_ERK;
    x_rdata[91] = B_Raf_Rap1_GTP_pMEK_ERK;
    x_rdata[92] = ppERK_MKP3;
    x_rdata[93] = dppERK_MKP3;
    x_rdata[94] = pro_TrkA;
    x_rdata[95] = NGF;
    x_rdata[96] = EGF;
    x_rdata[97] = pro_EGFR;
    x_rdata[98] = degradation;
}