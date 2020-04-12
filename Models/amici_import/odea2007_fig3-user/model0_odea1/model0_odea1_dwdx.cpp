#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_odea1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*IKK_deg_k_IKK_deg;
    dwdx[1] = 1.0*IkBaNFkB_cytoplasm*int_2ani_a7;
    dwdx[2] = 1.0*IkBbNFkB_cytoplasm*int_2bni_a8;
    dwdx[3] = 1.0*IkBeNFkB_cytoplasm*int_2eni_a9;
    dwdx[4] = 1.0*IkBa_cytoplasm*int_ai_a1;
    dwdx[5] = 1.0*IkBb_cytoplasm*int_bi_a2;
    dwdx[6] = 1.0*IkBe_cytoplasm*int_ei_a3;
    dwdx[7] = 1.0*deg_ai_r1;
    dwdx[8] = 1.0*NFkB_cytoplasm*int_2ain_a4_3;
    dwdx[9] = -1.0*int_ai_d1_1;
    dwdx[10] = 1.0*deg_ain_r4;
    dwdx[11] = -1.0*int_2ain_d4_3;
    dwdx[12] = -1.0*int_2ani_d1_2;
    dwdx[13] = 1.0*deg_an_deg4_c;
    dwdx[14] = 1.0*IKK*int_2ani_a7;
    dwdx[15] = -1.0*int_an_d4_1;
    dwdx[16] = 1.0*deg_an_n_deg4_n;
    dwdx[17] = -1.0*int_an_n_d4_2;
    dwdx[18] = 1.0*loc_an_k2_a;
    dwdx[19] = 1.0*deg_a_deg1_c;
    dwdx[20] = 1.0*IKK*int_ai_a1;
    dwdx[21] = 1.0*NFkB_cytoplasm*int_an_a4_1;
    dwdx[22] = 1.0*loc_a_tp1a;
    dwdx[23] = 1.0*mdeg_a_tr3a;
    dwdx[24] = 1.0*tsl_a_tr1a;
    dwdx[25] = 1.0*deg_a_n_deg1_n;
    dwdx[26] = 1.0*NFkB_nucleus*int_an_n_a4_2;
    dwdx[27] = -1.0*loc_a_tp2a;
    dwdx[28] = 1.0*deg_bi_r2;
    dwdx[29] = 1.0*NFkB_cytoplasm*int_2bin_a5_3;
    dwdx[30] = -1.0*int_bi_d2_1;
    dwdx[31] = 1.0*deg_bin_r5;
    dwdx[32] = -1.0*int_2bin_d5_3;
    dwdx[33] = -1.0*int_2bni_d2_2;
    dwdx[34] = 1.0*deg_bn_deg5_c;
    dwdx[35] = 1.0*IKK*int_2bni_a8;
    dwdx[36] = -1.0*int_bn_d5_1;
    dwdx[37] = 1.0*deg_bn_n_deg5_n;
    dwdx[38] = -1.0*int_bn_n_d5_2;
    dwdx[39] = 1.0*loc_bn_k2_b;
    dwdx[40] = 1.0*deg_b_deg2_c;
    dwdx[41] = 1.0*IKK*int_bi_a2;
    dwdx[42] = 1.0*NFkB_cytoplasm*int_bn_a5_1;
    dwdx[43] = 1.0*loc_b_tp1b;
    dwdx[44] = 1.0*mdeg_b_tr3b;
    dwdx[45] = 1.0*tsl_b_tr1b;
    dwdx[46] = 1.0*deg_b_n_deg2_n;
    dwdx[47] = 1.0*NFkB_nucleus*int_bn_n_a5_2;
    dwdx[48] = -1.0*loc_b_tp2b;
    dwdx[49] = 1.0*deg_ei_r3;
    dwdx[50] = 1.0*NFkB_cytoplasm*int_2ein_a6_3;
    dwdx[51] = -1.0*int_ei_d3_1;
    dwdx[52] = 1.0*deg_ein_r6;
    dwdx[53] = -1.0*int_2ein_d6_3;
    dwdx[54] = -1.0*int_2eni_d3_2;
    dwdx[55] = 1.0*deg_en_deg6_c;
    dwdx[56] = 1.0*IKK*int_2eni_a9;
    dwdx[57] = -1.0*int_en_d6_1;
    dwdx[58] = 1.0*deg_en_n_deg6_n;
    dwdx[59] = -1.0*int_en_n_d6_2;
    dwdx[60] = 1.0*loc_en_k2_e;
    dwdx[61] = 1.0*deg_e_deg3_c;
    dwdx[62] = 1.0*IKK*int_ei_a3;
    dwdx[63] = 1.0*NFkB_cytoplasm*int_en_a6_1;
    dwdx[64] = 1.0*loc_e_tp1e;
    dwdx[65] = 1.0*mdeg_e_tr3e;
    dwdx[66] = 1.0*tsl_e_tr1e;
    dwdx[67] = 1.0*deg_e_n_deg3_n;
    dwdx[68] = 1.0*NFkB_nucleus*int_en_n_a6_2;
    dwdx[69] = -1.0*loc_e_tp2e;
    dwdx[70] = 1.0*IkBaIKK*int_2ain_a4_3;
    dwdx[71] = 1.0*IkBbIKK*int_2bin_a5_3;
    dwdx[72] = 1.0*IkBeIKK*int_2ein_a6_3;
    dwdx[73] = 1.0*IkBa_cytoplasm*int_an_a4_1;
    dwdx[74] = 1.0*IkBb_cytoplasm*int_bn_a5_1;
    dwdx[75] = 1.0*IkBe_cytoplasm*int_en_a6_1;
    dwdx[76] = 1.0*loc_n_k1_2;
    dwdx[77] = 1.0*IkBa_nucleus*int_an_n_a4_2;
    dwdx[78] = 1.0*IkBb_nucleus*int_bn_n_a5_2;
    dwdx[79] = 1.0*IkBe_nucleus*int_en_n_a6_2;
    dwdx[80] = 2.0*NFkB_nucleus*itxn_a_tr2a_i;
    dwdx[81] = -1.0*loc_n_k1_1;
}