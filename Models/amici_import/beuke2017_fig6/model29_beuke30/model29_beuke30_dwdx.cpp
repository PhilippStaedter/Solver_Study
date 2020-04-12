#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model29_beuke30(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.2999999999999999e-5*IkBa_p65_association_k*p65_cytoplasm;
    dwdx[1] = 1.2999999999999999e-5*IkBa_deg_cyt_k;
    dwdx[2] = IkBa_nuclear_import_k;
    dwdx[3] = -1.2999999999999999e-5*IkBa_cytoplasm*IkBa_p65_phosphorylation_kcat*s15/pow(IkBa_cytoplasm + IkBa_phosphorylation_kM, 2) + 1.2999999999999999e-5*IkBa_p65_phosphorylation_kcat*s15/(IkBa_cytoplasm + IkBa_phosphorylation_kM);
    dwdx[4] = 1.2999999999999999e-5*IkBa_mRNA_deg_cyt_k;
    dwdx[5] = -1.2999999999999999e-5*IkBa_mRNA_cytoplasm*IkBa_translation_vmax/pow(IkBa_mRNA_cytoplasm + IkBa_translation_Km, 2) + 1.2999999999999999e-5*IkBa_translation_vmax/(IkBa_mRNA_cytoplasm + IkBa_translation_Km);
    dwdx[6] = 7.9999999999999996e-7*IkBa_mRNA_deg_nuc_k;
    dwdx[7] = IkBa_mRNA_transport_k;
    dwdx[8] = 7.9999999999999996e-7*IkBa_transcription_elongation_kbasal;
    dwdx[9] = 7.9999999999999996e-7*IkBa_transcription_elongation_kbasal;
    dwdx[10] = -1.2999999999999999e-5*JNK*JNK_active_phosphorylation_kcat*TNFa_TNFR1_EL/pow(JNK + JNK_active_phoshorylation_Km, 2) + 1.2999999999999999e-5*JNK_active_phosphorylation_kcat*TNFa_TNFR1_EL/(JNK + JNK_active_phoshorylation_Km);
    dwdx[11] = 1.2999999999999999e-5*JNK_basal_phosphorylation_k;
    dwdx[12] = 1.2999999999999999e-5*JNK_dephosphorylation_k;
    dwdx[13] = 2.5900000000000002e-6*LPS_LSEC_degradation_k;
    dwdx[14] = 2.5900000000000002e-6*LPS_MC_degradation_k;
    dwdx[15] = -4.7e-7*LPS*TNFa_LSEC_transcription_initiation_vA_LPS/pow(LPS + TNFa_LSEC_transcription_initiation_kA_LPS, 2) + 4.7e-7*TNFa_LSEC_transcription_initiation_vA_LPS/(LPS + TNFa_LSEC_transcription_initiation_kA_LPS);
    dwdx[16] = -3.4999999999999998e-7*LPS*TNFa_MC_transcription_initiation_vA_LPS/pow(LPS + TNFa_MC_transcription_initiation_kA_LPS, 2) + 3.4999999999999998e-7*TNFa_MC_transcription_initiation_vA_LPS/(LPS + TNFa_MC_transcription_initiation_kA_LPS);
    dwdx[17] = -7.9999999999999996e-7*MSK1*MSK1_phosphorylation_kcat*p38_P/pow(MSK1 + MSK1_activation_kM, 2) + 7.9999999999999996e-7*MSK1_phosphorylation_kcat*p38_P/(MSK1 + MSK1_activation_kM);
    dwdx[18] = 7.9999999999999996e-7*MSK1_dephosphorylation_k;
    dwdx[19] = 7.9999999999999996e-7*p65_MSK1_phosphrylation_kcat*s8/(p65_MSK1_phosphorylation_kM + s8);
    dwdx[20] = 7.9999999999999996e-7*p65_MSK1_phosphrylation_kcat*p65_nucleus/(p65_MSK1_phosphorylation_kM + p65_nucleus);
    dwdx[21] = -k_TNFR1_outer_membrane2vessicle_shuttle;
    dwdx[22] = 2.5900000000000002e-6*k_TNFa_TNFR1_association*s16;
    dwdx[23] = k_TNFR1_vessicle2outer_membrane_shuttle;
    dwdx[24] = -1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*TNFa_TNFR1_EL*s14/(pow(IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL, 2)*(IKKb_active_phosphorylation_kM + s14)) + 1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*s14/((IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL)*(IKKb_active_phosphorylation_kM + s14));
    dwdx[25] = 1.2999999999999999e-5*JNK*JNK_active_phosphorylation_kcat/(JNK + JNK_active_phoshorylation_Km);
    dwdx[26] = -4.9210000000000004e-11*k_TNFa_TNFR1_association;
    dwdx[27] = k_TNFa_TNFR1_internalisation;
    dwdx[28] = 1.2999999999999999e-5*p38*p38_active_phosphorylation_kcat/(p38 + p38_active_phoshorylation_Km);
    dwdx[29] = 1.2999999999999999e-5*k_TNFa_internal_degradation;
    dwdx[30] = TNFa_LSEC_translation_k;
    dwdx[31] = 4.7e-7*TNFa_mRNA_LSEC_degradation_k;
    dwdx[32] = TNFa_MC_translation_k;
    dwdx[33] = 3.4999999999999998e-7*TNFa_mRNA_MC_degradation_k;
    dwdx[34] = 4.7e-7*TNFa_LSEC_transcription_elongation_k;
    dwdx[35] = 3.4999999999999998e-7*TNFa_MC_transcription_elongation_k;
    dwdx[36] = -1.2999999999999999e-5*TNFa_TNFR1_EL*p38*p38_active_phosphorylation_kcat/pow(p38 + p38_active_phoshorylation_Km, 2) + 1.2999999999999999e-5*TNFa_TNFR1_EL*p38_active_phosphorylation_kcat/(p38 + p38_active_phoshorylation_Km);
    dwdx[37] = 1.2999999999999999e-5*p38_basal_phosphorylation_k;
    dwdx[38] = 7.9999999999999996e-7*MSK1*MSK1_phosphorylation_kcat/(MSK1 + MSK1_activation_kM);
    dwdx[39] = 1.2999999999999999e-5*p38_dephosphorylation_k;
    dwdx[40] = 7.9999999999999996e-7*p65_degradation_k;
    dwdx[41] = 7.9999999999999996e-7*p65_P_dephosphorylation_k;
    dwdx[42] = 7.9999999999999996e-7*IkBa_transcription_kA3_p65_2P;
    dwdx[43] = 7.9999999999999996e-7*IkBa_deg_complex_nuc_k;
    dwdx[44] = 7.9999999999999996e-7*p65_degradation_k;
    dwdx[45] = IkBa_p65_nuclear_export_k;
    dwdx[46] = -2.4e-9*IkBa_p65_association_k;
    dwdx[47] = 1.2999999999999999e-5*IkBa_cytoplasm*IkBa_p65_association_k;
    dwdx[48] = p65_nuclear_import_k;
    dwdx[49] = 1.2999999999999999e-5*p65_translation_k;
    dwdx[50] = -0.02*p65_nuclear_import_k;
    dwdx[51] = -7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat*p65_nucleus/pow(p65_MSK1_phosphorylation_kM + p65_nucleus, 2) + 7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat/(p65_MSK1_phosphorylation_kM + p65_nucleus);
    dwdx[52] = 7.9999999999999996e-7*IkBa_transcription_kA_p65;
    dwdx[53] = 7.9999999999999996e-7*IkBa_p65_association_k*s11;
    dwdx[54] = 7.9999999999999996e-7*IkBa_deg_nuc_k;
    dwdx[55] = -0.5*IkBa_nuclear_import_k;
    dwdx[56] = 7.9999999999999996e-7*IkBa_p65_association_k*p65_nucleus;
    dwdx[57] = -1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*TNFa_TNFR1_EL*s14/((IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL)*pow(IKKb_active_phosphorylation_kM + s14, 2)) + 1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*TNFa_TNFR1_EL/((IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL)*(IKKb_active_phosphorylation_kM + s14));
    dwdx[58] = 1.2999999999999999e-5*IKKb_basal_phosphoylation_k;
    dwdx[59] = 1.2999999999999999e-5*IkBa_cytoplasm*IkBa_p65_phosphorylation_kcat/(IkBa_cytoplasm + IkBa_phosphorylation_kM);
    dwdx[60] = 1.2999999999999999e-5*IkBa_p65_phosphorylation_kcat*s19/(IkBa_p65_phosphorylation_kM + s19);
    dwdx[61] = 1.2999999999999999e-5*IKKb_dephopshorylation_k;
    dwdx[62] = 2.5900000000000002e-6*TNFR1_EL*k_TNFa_TNFR1_association;
    dwdx[63] = 2.5900000000000002e-6*TNFa_degradation_k;
    dwdx[64] = -3.8999999999999998e-8*IkBa_p65_association_k;
    dwdx[65] = 1.2999999999999999e-5*IkBa_deg_complex_cyt_k;
    dwdx[66] = 1.2999999999999999e-5*p65_degradation_k;
    dwdx[67] = -IkBa_p65_nuclear_import_k;
    dwdx[68] = -1.2999999999999999e-5*IkBa_p65_phosphorylation_kcat*s15*s19/pow(IkBa_p65_phosphorylation_kM + s19, 2) + 1.2999999999999999e-5*IkBa_p65_phosphorylation_kcat*s15/(IkBa_p65_phosphorylation_kM + s19);
    dwdx[69] = 1.2999999999999999e-5*p65_degradation_k;
    dwdx[70] = 1.2999999999999999e-5*p65_P_dephosphorylation_k;
    dwdx[71] = p65_nuclear_import_k;
    dwdx[72] = 1.2999999999999999e-5*IkBa_p_active_degradation_k;
    dwdx[73] = 7.9999999999999996e-7*p65_degradation_k;
    dwdx[74] = 7.9999999999999996e-7*p65_P_dephosphorylation_k;
    dwdx[75] = -7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat*s8/pow(p65_MSK1_phosphorylation_kM + s8, 2) + 7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat/(p65_MSK1_phosphorylation_kM + s8);
    dwdx[76] = 7.9999999999999996e-7*IkBa_transcription_kA2_p65_P;
    dwdx[77] = -0.02*p65_nuclear_import_k;
}