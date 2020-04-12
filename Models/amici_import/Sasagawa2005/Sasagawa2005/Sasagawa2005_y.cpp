#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Sasagawa2005(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EGFR;
    y[1] = L_EGFR;
    y[2] = L_EGFR_dimer;
    y[3] = SOS;
    y[4] = L_dpEGFR;
    y[5] = pSOS;
    y[6] = SOS_Grb2;
    y[7] = Grb2;
    y[8] = Dok;
    y[9] = pDok;
    y[10] = Crk;
    y[11] = FRS2;
    y[12] = Shc;
    y[13] = pSOS_Grb2;
    y[14] = Rap1_GDP;
    y[15] = MEK;
    y[16] = MKP3;
    y[17] = pShc_dpEGFR;
    y[18] = dpEGFR_c_Cbl;
    y[19] = B_Raf_Rap1_GTP;
    y[20] = pShc_dpEGFR_c_Cbl;
    y[21] = pFRS2_dpEGFR_c_Cbl;
    y[22] = Shc_dpEGFR;
    y[23] = c_Cbl;
    y[24] = RasGAP;
    y[25] = c_Raf;
    y[26] = B_Raf;
    y[27] = ERK;
    y[28] = PP2A;
    y[29] = Ras_GDP;
    y[30] = Rap1GAP;
    y[31] = C3G;
    y[32] = NGFR;
    y[33] = pShc;
    y[34] = pFRS2_dpEGFR;
    y[35] = pTrkA_endo;
    y[36] = MEK_ERK;
    y[37] = pMEK_ERK;
    y[38] = FRS2_dpEGFR_c_Cbl_ubiq;
    y[39] = Crk_C3G_pFRS2_dpEGFR_c_Cbl;
    y[40] = pShc_dpEGFR_c_Cbl_ubiq;
    y[41] = Crk_C3G_pFRS2_dpEGFR;
    y[42] = Grb2_SOS_pShc_dpEGFR_c_Cbl_ubiq;
    y[43] = Grb2_SOS_pShc_dpEGFR_c_Cbl;
    y[44] = Shc_dpEGFR_c_Cbl_ubiq;
    y[45] = dpEGFR_c_Cbl_ubiq;
    y[46] = proteosome;
    y[47] = Grb2_SOS_pShc;
    y[48] = Shc_dpEGFR_c_Cbl;
    y[49] = Grb2_SOS_pShc_dpEGFR;
    y[50] = pFRS2;
    y[51] = FRS2_dpEGFR;
    y[52] = pDok_RasGAP;
    y[53] = pMEK;
    y[54] = FRS2_dpEGFR_c_Cbl;
    y[55] = pFRS2_dpEGFR_c_Cbl_ubiq;
    y[56] = Ras_GTP;
    y[57] = Crk_C3G_pFRS2_dpEGFR_c_Cbl_ubiq;
    y[58] = c_Raf_Ras_GTP;
    y[59] = B_Raf_Ras_GTP;
    y[60] = ppMEK;
    y[61] = ppERK;
    y[62] = pTrkA;
    y[63] = Crk_C3G;
    y[64] = Rap1_GTP;
    y[65] = L_NGFR;
    y[66] = ppMEK_ERK;
    y[67] = dppERK;
    y[68] = Shc_pTrkA;
    y[69] = Shc_pTrkA_endo;
    y[70] = pShc_pTrkA;
    y[71] = pFRS2_pTrkA;
    y[72] = FRS2_pTrkA;
    y[73] = pShc_pTrkA_endo;
    y[74] = FRS2_pTrkA_endo;
    y[75] = pFRS2_pTrkA_endo;
    y[76] = Crk_C3G_pFRS2_pTrkA_endo;
    y[77] = Grb2_SOS_pShc_pTrkA;
    y[78] = Crk_C3G_pFRS2_pTrkA;
    y[79] = Grb2_SOS_pShc_pTrkA_endo;
    y[80] = c_Raf_Ras_GTP_MEK;
    y[81] = c_Raf_Ras_GTP_pMEK;
    y[82] = c_Raf_Ras_GTP_MEK_ERK;
    y[83] = c_Raf_Ras_GTP_pMEK_ERK;
    y[84] = B_Raf_Ras_GTP_MEK;
    y[85] = B_Raf_Ras_GTP_pMEK;
    y[86] = B_Raf_Ras_GTP_MEK_ERK;
    y[87] = B_Raf_Ras_GTP_pMEK_ERK;
    y[88] = B_Raf_Rap1_GTP_MEK;
    y[89] = B_Raf_Rap1_GTP_pMEK;
    y[90] = B_Raf_Rap1_GTP_MEK_ERK;
    y[91] = B_Raf_Rap1_GTP_pMEK_ERK;
    y[92] = ppERK_MKP3;
    y[93] = dppERK_MKP3;
    y[94] = pro_TrkA;
    y[95] = NGF;
    y[96] = EGF;
    y[97] = pro_EGFR;
    y[98] = degradation;
}