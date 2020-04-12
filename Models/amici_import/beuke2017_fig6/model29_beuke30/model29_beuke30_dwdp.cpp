#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model29_beuke30(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[6] = -1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*TNFa_TNFR1_EL*s14/(pow(IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL, 2)*(IKKb_active_phosphorylation_kM + s14));
            break;
        case 1:
            dwdp[6] = -1.2999999999999999e-5*IKKb_active_phosphorylation_vmax*TNFa_TNFR1_EL*s14/((IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL)*pow(IKKb_active_phosphorylation_kM + s14, 2));
            break;
        case 2:
            dwdp[6] = 1.2999999999999999e-5*TNFa_TNFR1_EL*s14/((IKKb_active_phosphorylation_kA + TNFa_TNFR1_EL)*(IKKb_active_phosphorylation_kM + s14));
            break;
        case 3:
            dwdp[37] = 1.2999999999999999e-5*s14;
            break;
        case 4:
            dwdp[54] = 1.2999999999999999e-5*s15;
            break;
        case 5:
            dwdp[1] = 1.2999999999999999e-5*s19;
            break;
        case 6:
            dwdp[2] = 7.9999999999999996e-7*p65_IkBa_nucleus;
            break;
        case 7:
            dwdp[7] = 1.2999999999999999e-5*IkBa_cytoplasm;
            break;
        case 8:
            dwdp[8] = 7.9999999999999996e-7*s11;
            break;
        case 9:
            dwdp[9] = 1.2999999999999999e-5*IkBa_mRNA_cytoplasm;
            break;
        case 10:
            dwdp[10] = 7.9999999999999996e-7*IkBa_mRNA_nucleus;
            break;
        case 11:
            dwdp[11] = IkBa_mRNA_nucleus;
            break;
        case 12:
            dwdp[12] = IkBa_cytoplasm - 0.5*s11;
            break;
        case 13:
            dwdp[0] = 1.2999999999999999e-5*IkBa_cytoplasm*p65_cytoplasm - 3.8999999999999998e-8*s19;
            dwdp[55] = -2.4e-9*p65_IkBa_nucleus + 7.9999999999999996e-7*p65_nucleus*s11;
            break;
        case 14:
            dwdp[5] = p65_IkBa_nucleus;
            break;
        case 15:
            dwdp[5] = -s19;
            break;
        case 16:
            dwdp[53] = -1.2999999999999999e-5*IkBa_p65_phosphorylation_kcat*s15*s19/pow(IkBa_p65_phosphorylation_kM + s19, 2);
            break;
        case 17:
            dwdp[13] = 1.2999999999999999e-5*IkBa_cytoplasm*s15/(IkBa_cytoplasm + IkBa_phosphorylation_kM);
            dwdp[53] = 1.2999999999999999e-5*s15*s19/(IkBa_p65_phosphorylation_kM + s19);
            break;
        case 18:
            dwdp[52] = 1.2999999999999999e-5*s5;
            break;
        case 19:
            dwdp[13] = -1.2999999999999999e-5*IkBa_cytoplasm*IkBa_p65_phosphorylation_kcat*s15/pow(IkBa_cytoplasm + IkBa_phosphorylation_kM, 2);
            break;
        case 20:
            dwdp[14] = 7.9999999999999996e-7*IkBa_pre_mRNA_1;
            dwdp[15] = 7.9999999999999996e-7*IkBa_pre_mRNA_2;
            break;
        case 21:
            dwdp[51] = 7.9999999999999996e-7;
            break;
        case 22:
            dwdp[51] = 7.9999999999999996e-7*s8;
            break;
        case 23:
            dwdp[51] = 7.9999999999999996e-7*p65_2P;
            break;
        case 24:
            dwdp[51] = 7.9999999999999996e-7*p65_nucleus;
            break;
        case 25:
            dwdp[16] = -1.2999999999999999e-5*IkBa_mRNA_cytoplasm*IkBa_translation_vmax/pow(IkBa_mRNA_cytoplasm + IkBa_translation_Km, 2);
            break;
        case 26:
            dwdp[16] = 1.2999999999999999e-5*IkBa_mRNA_cytoplasm/(IkBa_mRNA_cytoplasm + IkBa_translation_Km);
            break;
        case 27:
            dwdp[17] = -1.2999999999999999e-5*JNK*JNK_active_phosphorylation_kcat*TNFa_TNFR1_EL/pow(JNK + JNK_active_phoshorylation_Km, 2);
            break;
        case 28:
            dwdp[17] = 1.2999999999999999e-5*JNK*TNFa_TNFR1_EL/(JNK + JNK_active_phoshorylation_Km);
            break;
        case 29:
            dwdp[18] = 1.2999999999999999e-5*JNK;
            break;
        case 30:
            dwdp[19] = 1.2999999999999999e-5*JNK_P;
            break;
        case 31:
            dwdp[20] = 2.5900000000000002e-6*LPS;
            break;
        case 32:
            dwdp[21] = 2.5900000000000002e-6*LPS;
            break;
        case 33:
            dwdp[22] = -7.9999999999999996e-7*MSK1*MSK1_phosphorylation_kcat*p38_P/pow(MSK1 + MSK1_activation_kM, 2);
            break;
        case 34:
            dwdp[23] = 7.9999999999999996e-7*MSK1_P;
            break;
        case 35:
            dwdp[22] = 7.9999999999999996e-7*MSK1*p38_P/(MSK1 + MSK1_activation_kM);
            break;
        case 36:
            dwdp[25] = 4.7e-7*TNFa_pre_mRNA_LSEC;
            break;
        case 37:
            dwdp[26] = -4.7e-7*LPS*TNFa_LSEC_transcription_initiation_vA_LPS/pow(LPS + TNFa_LSEC_transcription_initiation_kA_LPS, 2);
            break;
        case 38:
            dwdp[26] = 4.7e-7;
            break;
        case 39:
            dwdp[26] = 4.7e-7*LPS/(LPS + TNFa_LSEC_transcription_initiation_kA_LPS);
            break;
        case 40:
            dwdp[27] = TNFa_mRNA_LSEC;
            break;
        case 41:
            dwdp[28] = 3.4999999999999998e-7*TNFa_pre_mRNA_MC;
            break;
        case 42:
            dwdp[29] = -3.4999999999999998e-7*LPS*TNFa_MC_transcription_initiation_vA_LPS/pow(LPS + TNFa_MC_transcription_initiation_kA_LPS, 2);
            break;
        case 43:
            dwdp[29] = 3.4999999999999998e-7;
            break;
        case 44:
            dwdp[29] = 3.4999999999999998e-7*LPS/(LPS + TNFa_MC_transcription_initiation_kA_LPS);
            break;
        case 45:
            dwdp[30] = TNFa_mRNA_MC;
            break;
        case 46:
            dwdp[33] = 2.5900000000000002e-6*s16;
            break;
        case 47:
            dwdp[35] = 4.7e-7*TNFa_mRNA_LSEC;
            break;
        case 48:
            dwdp[36] = 3.4999999999999998e-7*TNFa_mRNA_MC;
            break;
        case 54:
            dwdp[24] = -TNFR1_EL;
            break;
        case 55:
            dwdp[24] = TNFR1_cytoplasm;
            break;
        case 56:
            dwdp[31] = 2.5900000000000002e-6*TNFR1_EL*s16 - 4.9210000000000004e-11*TNFa_TNFR1_EL;
            break;
        case 57:
            dwdp[32] = TNFa_TNFR1_EL;
            break;
        case 58:
            dwdp[34] = 1.2999999999999999e-5*TNFa_TNFR1_cytoplasm;
            break;
        case 62:
            dwdp[38] = -1.2999999999999999e-5*TNFa_TNFR1_EL*p38*p38_active_phosphorylation_kcat/pow(p38 + p38_active_phoshorylation_Km, 2);
            break;
        case 63:
            dwdp[38] = 1.2999999999999999e-5*TNFa_TNFR1_EL*p38/(p38 + p38_active_phoshorylation_Km);
            break;
        case 64:
            dwdp[39] = 1.2999999999999999e-5*p38;
            break;
        case 65:
            dwdp[40] = 1.2999999999999999e-5*p38_P;
            break;
        case 66:
            dwdp[47] = -7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat*s8/pow(p65_MSK1_phosphorylation_kM + s8, 2);
            dwdp[49] = -7.9999999999999996e-7*MSK1_P*p65_MSK1_phosphrylation_kcat*p65_nucleus/pow(p65_MSK1_phosphorylation_kM + p65_nucleus, 2);
            break;
        case 67:
            dwdp[47] = 7.9999999999999996e-7*MSK1_P*s8/(p65_MSK1_phosphorylation_kM + s8);
            dwdp[49] = 7.9999999999999996e-7*MSK1_P*p65_nucleus/(p65_MSK1_phosphorylation_kM + p65_nucleus);
            break;
        case 68:
            dwdp[42] = 7.9999999999999996e-7*p65_2P;
            dwdp[45] = 7.9999999999999996e-7*s8;
            dwdp[46] = 1.2999999999999999e-5*s3;
            break;
        case 69:
            dwdp[3] = 1.2999999999999999e-5*s19;
            dwdp[4] = 7.9999999999999996e-7*p65_IkBa_nucleus;
            dwdp[41] = 7.9999999999999996e-7*p65_2P;
            dwdp[43] = 1.2999999999999999e-5*s3;
            dwdp[44] = 7.9999999999999996e-7*s8;
            break;
        case 70:
            dwdp[48] = p65_cytoplasm - 0.02*p65_nucleus;
            dwdp[56] = s3 - 0.02*s8;
            break;
        case 71:
            dwdp[50] = 1.2999999999999999e-5*p65_mRNA;
            break;
    }
}