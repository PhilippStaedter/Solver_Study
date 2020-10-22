#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Pritchard2002(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = HXT_Vmax_1*(-GLCi + GLCo)/(HXT_Kglc_1*(GLCi*GLCo*HXT_Ki_1/pow(HXT_Kglc_1, 2) + 1 + (GLCi + GLCo)/HXT_Kglc_1));
    w[1] = 1.0*HK_Vmax_2*(-ADP*G6P/(HK_Katp_2*HK_Keq_2*HK_Kglc_2) + ATP*GLCi/(HK_Katp_2*HK_Kglc_2))/((ADP/HK_Kadp_2 + ATP/HK_Katp_2 + 1)*(G6P/HK_Kg6p_2 + GLCi/HK_Kglc_2 + 1));
    w[2] = 1.0*PGI_Vmax_3*(-F6P/(PGI_Keq_3*PGI_Kg6p_3) + G6P/PGI_Kg6p_3)/(F6P/PGI_Kf6p_3 + G6P/PGI_Kg6p_3 + 1);
    w[3] = 1.0*ATP*F6P*PFK_Vmax_4*PFK_gR_4*(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1)/(PFK_Katp_4*PFK_Kf6p_4*(PFK_L0_4*pow(AMP*PFK_Camp_4/PFK_Kamp_4 + 1, 2)*pow(ATP*PFK_Catp_4/PFK_Katp_4 + 1, 2)*pow(ATP*PFK_Ciatp_4/PFK_Kiatp_4 + 1, 2)*pow(F16bP*PFK_Cf16_4/PFK_Kf16_4 + F26bP*PFK_Cf26_4/PFK_Kf26_4 + 1, 2)/(pow(AMP/PFK_Kamp_4 + 1, 2)*pow(ATP/PFK_Kiatp_4 + 1, 2)*pow(F16bP/PFK_Kf16_4 + F26bP/PFK_Kf26_4 + 1, 2)) + pow(ATP*F6P*PFK_gR_4/(PFK_Katp_4*PFK_Kf6p_4) + ATP/PFK_Katp_4 + F6P/PFK_Kf6p_4 + 1, 2)));
    w[4] = 1.0*ALD_Vmax_5*(F16bP/ALD_Kf16bp_5 - DHAP*GAP/(ALD_Keq_5*ALD_Kf16bp_5))/(1 + GAP/ALD_Kgap_5 + F16bP/ALD_Kf16bp_5 + F16bP*GAP/(ALD_Kf16bp_5*ALD_Kigap_5) + DHAP/ALD_Kdhap_5 + DHAP*GAP/(ALD_Kdhap_5*ALD_Kgap_5));
    w[5] = 1.0*DHAP*TPI_k1_6 - 1.0*GAP*TPI_k2_6;
    w[6] = 1.0*GAPDH_C_7*(-BPG*GAPDH_Vmaxr_7*NADH/(GAPDH_Kbpg_7*GAPDH_Knadh_7) + GAP*GAPDH_Vmaxf_7*NAD/(GAPDH_Kgap_7*GAPDH_Knad_7))/((1 + NADH/GAPDH_Knadh_7 + NAD/GAPDH_Knad_7)*(BPG/GAPDH_Kbpg_7 + GAP/GAPDH_Kgap_7 + 1));
    w[7] = 1.0*PGK_Vmax_8*(ADP*BPG*PGK_Keq_8 - ATP*P3G)/(PGK_Katp_8*PGK_Kp3g_8*(ADP/PGK_Kadp_8 + ATP/PGK_Katp_8 + 1)*(BPG/PGK_Kbpg_8 + P3G/PGK_Kp3g_8 + 1));
    w[8] = 1.0*PGM_Vmax_9*(-P2G/(PGM_Keq_9*PGM_Kp3g_9) + P3G/PGM_Kp3g_9)/(P2G/PGM_Kp2g_9 + P3G/PGM_Kp3g_9 + 1);
    w[9] = 1.0*ENO_Vmax_10*(P2G/ENO_Kp2g_10 - PEP/(ENO_Keq_10*ENO_Kp2g_10))/(1 + PEP/ENO_Kpep_10 + P2G/ENO_Kp2g_10);
    w[10] = 1.0*PYK_Vmax_11*(ADP*PEP/(PYK_Kadp_11*PYK_Kpep_11) - ATP*PYR/(PYK_Kadp_11*PYK_Keq_11*PYK_Kpep_11))/((ADP/PYK_Kadp_11 + ATP/PYK_Katp_11 + 1)*(PEP/PYK_Kpep_11 + 1 + PYR/PYK_Kpyr_11));
    w[11] = 1.0*PDC_Vmax_12*pow(PYR/PDC_Kpyr_12, PDC_nH_12)/(pow(PYR/PDC_Kpyr_12, PDC_nH_12) + 1);
    w[12] = 1.0*ADH_Vmax_13*(EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) - AcAld*NADH/(ADH_Keq_13*ADH_Ketoh_13*ADH_Kinad_13))/(1 + NADH/ADH_Kinadh_13 + NAD/ADH_Kinad_13 + ADH_Knad_13*EtOH/(ADH_Ketoh_13*ADH_Kinad_13) + EtOH*NAD/(ADH_Ketoh_13*ADH_Kinad_13) + ADH_Knad_13*EtOH*NADH/(ADH_Ketoh_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NAD/(ADH_Ketoh_13*ADH_Kiacald_13*ADH_Kinad_13) + ADH_Knadh_13*AcAld/(ADH_Kacald_13*ADH_Kinadh_13) + AcAld*NADH/(ADH_Kacald_13*ADH_Kinadh_13) + ADH_Knadh_13*AcAld*NAD/(ADH_Kacald_13*ADH_Kinad_13*ADH_Kinadh_13) + AcAld*EtOH*NADH/(ADH_Kacald_13*ADH_Kietoh_13*ADH_Kinadh_13));
    w[13] = 1.0*ATP*ATPase_Katpase_14;
    w[14] = 1.0*pow(ADP, 2)*AK_k1_15 - 1.0*AK_k2_15*AMP*ATP;
    w[15] = 1.0*G3PDH_Vmax_16*(DHAP*NADH/(G3PDH_Kdhap_16*G3PDH_Knadh_16) - Glycerol*NAD/(G3PDH_Kdhap_16*G3PDH_Keq_16*G3PDH_Knadh_16))/((1 + NADH/G3PDH_Knadh_16 + NAD/G3PDH_Knad_16)*(DHAP/G3PDH_Kdhap_16 + 1 + Glycerol/G3PDH_Kglycerol_16));
    w[16] = 1.0*Glycogen_Branch_KGLYCOGEN_17;
    w[17] = 1.0*Trehalose_Branch_Ktrehalose_18;
    w[18] = 1.0*AcAld*Succinate_Branch_k_19;
}