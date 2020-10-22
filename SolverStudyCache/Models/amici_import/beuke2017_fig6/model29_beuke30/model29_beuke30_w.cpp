#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model29_beuke30(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.2999999999999999e-5*IkBa_cytoplasm*IkBa_p65_association_k*p65_cytoplasm - 3.8999999999999998e-8*IkBa_p65_association_k*s19;
    w[1] = 1.2999999999999999e-5*IkBa_deg_complex_cyt_k*s19;
    w[2] = 7.9999999999999996e-7*IkBa_deg_complex_nuc_k*p65_IkBa_nucleus;
    w[3] = 1.2999999999999999e-5*p65_degradation_k*s19;
    w[4] = 7.9999999999999996e-7*p65_IkBa_nucleus*p65_degradation_k;
    w[5] = IkBa_p65_nuclear_export_k*p65_IkBa_nucleus - IkBa_p65_nuclear_import_k*s19;
    w[6] = 1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*TNFa_TNFR1_EL*s14/((IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL)*(IKKb_active_phosphorylation_kM + s14));
    w[7] = 1.2999999999999999e-5*IkBa_cytoplasm*IkBa_deg_cyt_k;
    w[8] = 7.9999999999999996e-7*IkBa_deg_nuc_k*s11;
    w[9] = 1.2999999999999999e-5*IkBa_mRNA_cytoplasm*IkBa_mRNA_deg_cyt_k;
    w[10] = 7.9999999999999996e-7*IkBa_mRNA_deg_nuc_k*IkBa_mRNA_nucleus;
    w[11] = IkBa_mRNA_nucleus*IkBa_mRNA_transport_k;
    w[12] = IkBa_cytoplasm*IkBa_nuclear_import_k - 0.5*IkBa_nuclear_import_k*s11;
    w[13] = 1.2999999999999999e-5*IkBa_cytoplasm*IkBa_p65_phosphorylation_kcat*s15/(IkBa_cytoplasm + IkBa_phosphorylation_kM);
    w[14] = 7.9999999999999996e-7*IkBa_pre_mRNA_1*IkBa_transcription_elongation_kbasal;
    w[15] = 7.9999999999999996e-7*IkBa_pre_mRNA_2*IkBa_transcription_elongation_kbasal;
    w[16] = 1.2999999999999999e-5*IkBa_mRNA_cytoplasm*IkBa_translation_vmax/(IkBa_mRNA_cytoplasm + IkBa_translation_Km);
    w[17] = 1.2999999999999999e-5*JNK*JNK_active_phosphorylation_kcat*TNFa_TNFR1_EL/(JNK + JNK_active_phoshorylation_Km);
    w[18] = 1.2999999999999999e-5*JNK*JNK_basal_phosphorylation_k;
    w[19] = 1.2999999999999999e-5*JNK_P*JNK_dephosphorylation_k;
    w[20] = 2.5900000000000002e-6*LPS*LPS_LSEC_degradation_k;
    w[21] = 2.5900000000000002e-6*LPS*LPS_MC_degradation_k;
    w[22] = 7.9999999999999996e-7*MSK1*MSK1_phosphorylation_kcat*p38_P/(MSK1 + MSK1_activation_kM);
    w[23] = 7.9999999999999996e-7*MSK1_P*MSK1_dephosphorylation_k;
    w[24] = -TNFR1_EL*k_TNFR1_outer_membrane2vessicle_shuttle + TNFR1_cytoplasm*k_TNFR1_vessicle2outer_membrane_shuttle;
    w[25] = 4.7e-7*TNFa_LSEC_transcription_elongation_k*TNFa_pre_mRNA_LSEC;
    w[26] = 4.7e-7*LPS*TNFa_LSEC_transcription_initiation_vA_LPS/(LPS + TNFa_LSEC_transcription_initiation_kA_LPS) + 4.7e-7*TNFa_LSEC_transcription_initiation_kbasal;
    w[27] = TNFa_LSEC_translation_k*TNFa_mRNA_LSEC;
    w[28] = 3.4999999999999998e-7*TNFa_MC_transcription_elongation_k*TNFa_pre_mRNA_MC;
    w[29] = 3.4999999999999998e-7*LPS*TNFa_MC_transcription_initiation_vA_LPS/(LPS + TNFa_MC_transcription_initiation_kA_LPS) + 3.4999999999999998e-7*TNFa_MC_transcription_initiation_kbasal;
    w[30] = TNFa_MC_translation_k*TNFa_mRNA_MC;
    w[31] = 2.5900000000000002e-6*TNFR1_EL*k_TNFa_TNFR1_association*s16 - 4.9210000000000004e-11*TNFa_TNFR1_EL*k_TNFa_TNFR1_association;
    w[32] = TNFa_TNFR1_EL*k_TNFa_TNFR1_internalisation;
    w[33] = 2.5900000000000002e-6*TNFa_degradation_k*s16;
    w[34] = 1.2999999999999999e-5*TNFa_TNFR1_cytoplasm*k_TNFa_internal_degradation;
    w[35] = 4.7e-7*TNFa_mRNA_LSEC*TNFa_mRNA_LSEC_degradation_k;
    w[36] = 3.4999999999999998e-7*TNFa_mRNA_MC*TNFa_mRNA_MC_degradation_k;
    w[37] = 1.2999999999999999e-5*IKKb_basal_phosphoylation_k*s14;
    w[38] = 1.2999999999999999e-5*TNFa_TNFR1_EL*p38*p38_active_phosphorylation_kcat/(p38 + p38_active_phoshorylation_Km);
    w[39] = 1.2999999999999999e-5*p38*p38_basal_phosphorylation_k;
    w[40] = 1.2999999999999999e-5*p38_P*p38_dephosphorylation_k;
    w[41] = 7.9999999999999996e-7*p65_2P*p65_degradation_k;
    w[42] = 7.9999999999999996e-7*p65_2P*p65_P_dephosphorylation_k;
    w[43] = 1.2999999999999999e-5*p65_degradation_k*s3;
    w[44] = 7.9999999999999996e-7*p65_degradation_k*s8;
    w[45] = 7.9999999999999996e-7*p65_P_dephosphorylation_k*s8;
    w[46] = 1.2999999999999999e-5*p65_P_dephosphorylation_k*s3;
    w[47] = 7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat*s8/(p65_MSK1_phosphorylation_kM + s8);
    w[48] = p65_cytoplasm*p65_nuclear_import_k - 0.02*p65_nuclear_import_k*p65_nucleus;
    w[49] = 7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat*p65_nucleus/(p65_MSK1_phosphorylation_kM + p65_nucleus);
    w[50] = 1.2999999999999999e-5*p65_mRNA*p65_translation_k;
    w[51] = 7.9999999999999996e-7*IkBa_transcription_initiation_kbasal + 7.9999999999999996e-7*IkBa_transcription_kA2_p65_P*s8 + 7.9999999999999996e-7*IkBa_transcription_kA3_p65_2P*p65_2P + 7.9999999999999996e-7*IkBa_transcription_kA_p65*p65_nucleus;
    w[52] = 1.2999999999999999e-5*IkBa_p_active_degradation_k*s5;
    w[53] = 1.2999999999999999e-5*IkBa_p65_phosphorylation_kcat*s15*s19/(IkBa_p65_phosphorylation_kM + s19);
    w[54] = 1.2999999999999999e-5*IKKb_dephopshorylation_k*s15;
    w[55] = -2.4e-9*IkBa_p65_association_k*p65_IkBa_nucleus + 7.9999999999999996e-7*IkBa_p65_association_k*p65_nucleus*s11;
    w[56] = p65_nuclear_import_k*s3 - 0.02*p65_nuclear_import_k*s8;
}