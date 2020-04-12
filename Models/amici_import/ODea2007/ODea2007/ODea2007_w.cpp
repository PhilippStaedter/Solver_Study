#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_ODea2007(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*txn_a_tr2a;
    w[1] = 1.0*IkBa_mRNA*mdeg_a_tr3a;
    w[2] = 1.0*IkBa_mRNA*tsl_a_tr1a;
    w[3] = 1.0*IKK*IkBa_cytoplasm*int_ai_a1 - 1.0*IkBaIKK*int_ai_d1_1;
    w[4] = -1.0*IkBaNFkB_cytoplasm*int_an_d4_1 + 1.0*IkBa_cytoplasm*NFkB_cytoplasm*int_an_a4_1;
    w[5] = -1.0*IkBaNFkB_nucleus*int_an_n_d4_2 + 1.0*IkBa_nucleus*NFkB_nucleus*int_an_n_a4_2;
    w[6] = 1.0*IKK*IkBaNFkB_cytoplasm*int_2ani_a7 - 1.0*IkBaIKKNFkB*int_2ani_d1_2;
    w[7] = 1.0*IkBaIKK*NFkB_cytoplasm*int_2ain_a4_3 - 1.0*IkBaIKKNFkB*int_2ain_d4_3;
    w[8] = 1.0*IkBa_cytoplasm*deg_a_deg1_c;
    w[9] = 1.0*IkBa_nucleus*deg_a_n_deg1_n;
    w[10] = 1.0*IkBaNFkB_nucleus*deg_an_n_deg4_n;
    w[11] = 1.0*IkBaNFkB_cytoplasm*deg_an_deg4_c;
    w[12] = 1.0*IkBaIKK*deg_ai_r1;
    w[13] = 1.0*IkBaIKKNFkB*deg_ain_r4;
    w[14] = 1.0*IkBa_cytoplasm*loc_a_tp1a - 1.0*IkBa_nucleus*loc_a_tp2a;
    w[15] = 1.0*IkBaNFkB_nucleus*loc_an_k2_a;
    w[16] = 1.0*IkBbNFkB_nucleus*loc_bn_k2_b;
    w[17] = 1.0*IkBb_cytoplasm*loc_b_tp1b - 1.0*IkBb_nucleus*loc_b_tp2b;
    w[18] = 1.0*IkBbIKKNFkB*deg_bin_r5;
    w[19] = 1.0*IkBbIKK*deg_bi_r2;
    w[20] = 1.0*IkBbNFkB_cytoplasm*deg_bn_deg5_c;
    w[21] = 1.0*IkBbNFkB_nucleus*deg_bn_n_deg5_n;
    w[22] = 1.0*IkBb_nucleus*deg_b_n_deg2_n;
    w[23] = 1.0*IkBb_cytoplasm*deg_b_deg2_c;
    w[24] = 1.0*IkBbIKK*NFkB_cytoplasm*int_2bin_a5_3 - 1.0*IkBbIKKNFkB*int_2bin_d5_3;
    w[25] = 1.0*IKK*IkBbNFkB_cytoplasm*int_2bni_a8 - 1.0*IkBbIKKNFkB*int_2bni_d2_2;
    w[26] = -1.0*IkBbNFkB_nucleus*int_bn_n_d5_2 + 1.0*IkBb_nucleus*NFkB_nucleus*int_bn_n_a5_2;
    w[27] = -1.0*IkBbNFkB_cytoplasm*int_bn_d5_1 + 1.0*IkBb_cytoplasm*NFkB_cytoplasm*int_bn_a5_1;
    w[28] = 1.0*IKK*IkBb_cytoplasm*int_bi_a2 - 1.0*IkBbIKK*int_bi_d2_1;
    w[29] = 1.0*IkBb_mRNA*tsl_b_tr1b;
    w[30] = 1.0*IkBb_mRNA*mdeg_b_tr3b;
    w[31] = 1.0*txn_b_tr2b;
    w[32] = 1.0*NFkB_cytoplasm*loc_n_k1_2 - 1.0*NFkB_nucleus*loc_n_k1_1;
    w[33] = 1.0*txn_e_tr2e;
    w[34] = 1.0*IkBe_mRNA*mdeg_e_tr3e;
    w[35] = 1.0*IkBe_mRNA*tsl_e_tr1e;
    w[36] = 1.0*IKK*IkBe_cytoplasm*int_ei_a3 - 1.0*IkBeIKK*int_ei_d3_1;
    w[37] = -1.0*IkBeNFkB_cytoplasm*int_en_d6_1 + 1.0*IkBe_cytoplasm*NFkB_cytoplasm*int_en_a6_1;
    w[38] = -1.0*IkBeNFkB_nucleus*int_en_n_d6_2 + 1.0*IkBe_nucleus*NFkB_nucleus*int_en_n_a6_2;
    w[39] = 1.0*IKK*IkBeNFkB_cytoplasm*int_2eni_a9 - 1.0*IkBeIKKNFkB*int_2eni_d3_2;
    w[40] = 1.0*IkBeIKK*NFkB_cytoplasm*int_2ein_a6_3 - 1.0*IkBeIKKNFkB*int_2ein_d6_3;
    w[41] = 1.0*IkBe_cytoplasm*deg_e_deg3_c;
    w[42] = 1.0*IkBe_nucleus*deg_e_n_deg3_n;
    w[43] = 1.0*IkBeNFkB_nucleus*deg_en_n_deg6_n;
    w[44] = 1.0*IkBeNFkB_cytoplasm*deg_en_deg6_c;
    w[45] = 1.0*IkBeIKK*deg_ei_r3;
    w[46] = 1.0*IkBeIKKNFkB*deg_ein_r6;
    w[47] = 1.0*IkBe_cytoplasm*loc_e_tp1e - 1.0*IkBe_nucleus*loc_e_tp2e;
    w[48] = 1.0*IkBeNFkB_nucleus*loc_en_k2_e;
    w[49] = 1.0*IKK*IKK_deg_k_IKK_deg;
    w[50] = 1.0*pow(NFkB_nucleus, 2)*itxn_a_tr2a_i;
}