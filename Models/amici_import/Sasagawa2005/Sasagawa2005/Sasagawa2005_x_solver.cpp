#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Sasagawa2005(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EGFR;
    x_solver[1] = L_EGFR;
    x_solver[2] = L_EGFR_dimer;
    x_solver[3] = SOS;
    x_solver[4] = L_dpEGFR;
    x_solver[5] = pSOS;
    x_solver[6] = SOS_Grb2;
    x_solver[7] = Grb2;
    x_solver[8] = Dok;
    x_solver[9] = pDok;
    x_solver[10] = Crk;
    x_solver[11] = FRS2;
    x_solver[12] = Shc;
    x_solver[13] = pSOS_Grb2;
    x_solver[14] = Rap1_GDP;
    x_solver[15] = MEK;
    x_solver[16] = MKP3;
    x_solver[17] = pShc_dpEGFR;
    x_solver[18] = dpEGFR_c_Cbl;
    x_solver[19] = B_Raf_Rap1_GTP;
    x_solver[20] = pShc_dpEGFR_c_Cbl;
    x_solver[21] = pFRS2_dpEGFR_c_Cbl;
    x_solver[22] = Shc_dpEGFR;
    x_solver[23] = c_Cbl;
    x_solver[24] = RasGAP;
    x_solver[25] = c_Raf;
    x_solver[26] = B_Raf;
    x_solver[27] = ERK;
    x_solver[28] = PP2A;
    x_solver[29] = Ras_GDP;
    x_solver[30] = Rap1GAP;
    x_solver[31] = C3G;
    x_solver[32] = NGFR;
    x_solver[33] = pShc;
    x_solver[34] = pFRS2_dpEGFR;
    x_solver[35] = pTrkA_endo;
    x_solver[36] = MEK_ERK;
    x_solver[37] = pMEK_ERK;
    x_solver[38] = FRS2_dpEGFR_c_Cbl_ubiq;
    x_solver[39] = Crk_C3G_pFRS2_dpEGFR_c_Cbl;
    x_solver[40] = pShc_dpEGFR_c_Cbl_ubiq;
    x_solver[41] = Crk_C3G_pFRS2_dpEGFR;
    x_solver[42] = Grb2_SOS_pShc_dpEGFR_c_Cbl_ubiq;
    x_solver[43] = Grb2_SOS_pShc_dpEGFR_c_Cbl;
    x_solver[44] = Shc_dpEGFR_c_Cbl_ubiq;
    x_solver[45] = dpEGFR_c_Cbl_ubiq;
    x_solver[46] = proteosome;
    x_solver[47] = Grb2_SOS_pShc;
    x_solver[48] = Shc_dpEGFR_c_Cbl;
    x_solver[49] = Grb2_SOS_pShc_dpEGFR;
    x_solver[50] = pFRS2;
    x_solver[51] = FRS2_dpEGFR;
    x_solver[52] = pDok_RasGAP;
    x_solver[53] = pMEK;
    x_solver[54] = FRS2_dpEGFR_c_Cbl;
    x_solver[55] = pFRS2_dpEGFR_c_Cbl_ubiq;
    x_solver[56] = Ras_GTP;
    x_solver[57] = Crk_C3G_pFRS2_dpEGFR_c_Cbl_ubiq;
    x_solver[58] = c_Raf_Ras_GTP;
    x_solver[59] = B_Raf_Ras_GTP;
    x_solver[60] = ppMEK;
    x_solver[61] = ppERK;
    x_solver[62] = pTrkA;
    x_solver[63] = Crk_C3G;
    x_solver[64] = Rap1_GTP;
    x_solver[65] = L_NGFR;
    x_solver[66] = ppMEK_ERK;
    x_solver[67] = dppERK;
    x_solver[68] = Shc_pTrkA;
    x_solver[69] = Shc_pTrkA_endo;
    x_solver[70] = pShc_pTrkA;
    x_solver[71] = pFRS2_pTrkA;
    x_solver[72] = FRS2_pTrkA;
    x_solver[73] = pShc_pTrkA_endo;
    x_solver[74] = FRS2_pTrkA_endo;
    x_solver[75] = pFRS2_pTrkA_endo;
    x_solver[76] = Crk_C3G_pFRS2_pTrkA_endo;
    x_solver[77] = Grb2_SOS_pShc_pTrkA;
    x_solver[78] = Crk_C3G_pFRS2_pTrkA;
    x_solver[79] = Grb2_SOS_pShc_pTrkA_endo;
    x_solver[80] = c_Raf_Ras_GTP_MEK;
    x_solver[81] = c_Raf_Ras_GTP_pMEK;
    x_solver[82] = c_Raf_Ras_GTP_MEK_ERK;
    x_solver[83] = c_Raf_Ras_GTP_pMEK_ERK;
    x_solver[84] = B_Raf_Ras_GTP_MEK;
    x_solver[85] = B_Raf_Ras_GTP_pMEK;
    x_solver[86] = B_Raf_Ras_GTP_MEK_ERK;
    x_solver[87] = B_Raf_Ras_GTP_pMEK_ERK;
    x_solver[88] = B_Raf_Rap1_GTP_MEK;
    x_solver[89] = B_Raf_Rap1_GTP_pMEK;
    x_solver[90] = B_Raf_Rap1_GTP_MEK_ERK;
    x_solver[91] = B_Raf_Rap1_GTP_pMEK_ERK;
    x_solver[92] = ppERK_MKP3;
    x_solver[93] = dppERK_MKP3;
    x_solver[94] = pro_TrkA;
    x_solver[95] = NGF;
    x_solver[96] = EGF;
    x_solver[97] = pro_EGFR;
    x_solver[98] = degradation;
}