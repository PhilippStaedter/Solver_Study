#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model3_mcauley1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = multiplier*reaction_1_k1;
    dwdx[1] = reaction_31_k26*species_11*species_31;
    dwdx[2] = reaction_14_k11;
    dwdx[3] = -reaction_26_kprs*species_26/pow(species_11, 2);
    dwdx[4] = reaction_28_k23*species_14;
    dwdx[5] = reaction_30_k1;
    dwdx[6] = reaction_31_k26*species_10*species_31;
    dwdx[7] = -reaction_32_PCSS*reaction_32_PCSmax*pow(species_11/reaction_32_PPCt, reaction_32_PCSS)/(species_11*pow(pow(species_11/reaction_32_PPCt, reaction_32_PCSS) + 1, 2));
    dwdx[8] = reaction_9_k8;
    dwdx[9] = reaction_13_k10*species_15;
    dwdx[10] = reaction_12_k9*species_7;
    dwdx[11] = reaction_28_k23*species_11;
    dwdx[12] = reaction_13_k10*species_13;
    dwdx[13] = reaction_29_k24*species_28;
    dwdx[14] = reaction_18_k1;
    dwdx[15] = reaction_19_k15*species_22;
    dwdx[16] = reaction_17_k1;
    dwdx[17] = reaction_22_k18*species_23;
    dwdx[18] = reaction_16_khrs/species_7;
    dwdx[19] = -reaction_2_ICSmax*reaction_2_IS*pow(species_2/reaction_2_ICt, reaction_2_IS)/(species_2*pow(pow(species_2/reaction_2_ICt, reaction_2_IS) + 1, 2));
    dwdx[20] = reaction_7_k6*species_5;
    dwdx[21] = reaction_8_k7*species_5;
    dwdx[22] = reaction_20_k1;
    dwdx[23] = reaction_21_k17*species_24;
    dwdx[24] = reaction_19_k15*species_17;
    dwdx[25] = reaction_22_k18*species_18;
    dwdx[26] = reaction_23_k1;
    dwdx[27] = reaction_24_k20*species_25;
    dwdx[28] = reaction_25_k1;
    dwdx[29] = reaction_21_k17*species_21;
    dwdx[30] = reaction_24_k20*species_23;
    dwdx[31] = reaction_27_k1;
    dwdx[32] = reaction_26_kprs/species_11;
    dwdx[33] = reaction_29_k24*species_15;
    dwdx[34] = reaction_33_k27*species_33;
    dwdx[35] = reaction_34_k28*species_33;
    dwdx[36] = reaction_35_k29*species_34;
    dwdx[37] = reaction_31_k26*species_10*species_11;
    dwdx[38] = reaction_33_k27*species_30;
    dwdx[39] = reaction_34_k28*species_30;
    dwdx[40] = reaction_35_k29*species_30;
    dwdx[41] = reaction_3_k1;
    dwdx[42] = -reaction_6_k5*species_7/pow(species_4, 2);
    dwdx[43] = reaction_4_k1;
    dwdx[44] = reaction_5_k1;
    dwdx[45] = reaction_7_k6*species_2;
    dwdx[46] = reaction_8_k7*species_2;
    dwdx[47] = reaction_10_BCRmax*reaction_10_BS*pow(reaction_10_BCRt/species_7, reaction_10_BS)/(species_7*pow(pow(reaction_10_BCRt/species_7, reaction_10_BS) + 1, 2));
    dwdx[48] = -reaction_11_HCSmax*reaction_11_HS*pow(species_7/reaction_11_HCSt, reaction_11_HS)/(species_7*pow(pow(species_7/reaction_11_HCSt, reaction_11_HS) + 1, 2));
    dwdx[49] = reaction_12_k9*species_14;
    dwdx[50] = reaction_15_k1;
    dwdx[51] = -reaction_16_khrs*species_19/pow(species_7, 2);
    dwdx[52] = reaction_6_k5/species_4;
}