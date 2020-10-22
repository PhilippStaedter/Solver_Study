#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model3_mcauley1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = multiplier*reaction_1_k1*species_1;
    w[1] = reaction_10_BCRmax/(pow(reaction_10_BCRt/species_7, reaction_10_BS) + 1);
    w[2] = reaction_11_HCSmax/(pow(species_7/reaction_11_HCSt, reaction_11_HS) + 1);
    w[3] = reaction_12_k9*species_14*species_7;
    w[4] = reaction_13_k10*species_13*species_15;
    w[5] = reaction_14_k11*species_11;
    w[6] = reaction_15_k1*species_7;
    w[7] = reaction_16_khrs*species_19/species_7;
    w[8] = reaction_17_k1*species_18;
    w[9] = reaction_18_k1*species_17;
    w[10] = reaction_19_k15*species_17*species_22;
    w[11] = reaction_2_ICSmax/(pow(species_2/reaction_2_ICt, reaction_2_IS) + 1);
    w[12] = reaction_20_k1*species_21;
    w[13] = reaction_21_k17*species_21*species_24;
    w[14] = reaction_22_k18*species_18*species_23;
    w[15] = reaction_23_k1*species_23;
    w[16] = reaction_24_k20*species_23*species_25;
    w[17] = reaction_25_k1*species_23;
    w[18] = reaction_26_kprs*species_26/species_11;
    w[19] = reaction_27_k1*species_25;
    w[20] = reaction_28_k23*species_11*species_14;
    w[21] = reaction_29_k24*species_15*species_28;
    w[22] = reaction_3_k1*species_4;
    w[23] = reaction_30_k1*species_11;
    w[24] = reaction_31_k26*species_10*species_11*species_31;
    w[25] = reaction_32_PCSmax/(pow(species_11/reaction_32_PPCt, reaction_32_PCSS) + 1);
    w[26] = reaction_33_k27*species_30*species_33;
    w[27] = reaction_34_k28*species_30*species_33;
    w[28] = reaction_35_k29*species_30*species_34;
    w[29] = reaction_4_k1*species_5;
    w[30] = reaction_5_k1*species_5;
    w[31] = reaction_6_k5*species_7/species_4;
    w[32] = reaction_7_k6*species_2*species_5;
    w[33] = reaction_8_k7*species_2*species_5;
    w[34] = reaction_9_k8*species_11;
}