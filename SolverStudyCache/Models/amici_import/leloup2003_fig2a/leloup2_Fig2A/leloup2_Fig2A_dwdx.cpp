#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_leloup2_Fig2A(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*reaction_1_k;
    dwdx[1] = 1.0*reaction_2_k1;
    dwdx[2] = -1.0*reaction_24_V*species_0/pow(reaction_24_Km + species_0, 2) + 1.0*reaction_24_V/(reaction_24_Km + species_0);
    dwdx[3] = -1.0*reaction_3_V*species_1/pow(reaction_3_Km + species_1, 2) + 1.0*reaction_3_V/(reaction_3_Km + species_1);
    dwdx[4] = 1.0*reaction_34_k1;
    dwdx[5] = 1.0*reaction_4_k1;
    dwdx[6] = -1.0*reaction_13_k2;
    dwdx[7] = -1.0*reaction_14_V*species_10/pow(reaction_14_Km + species_10, 2) + 1.0*reaction_14_V/(reaction_14_Km + species_10);
    dwdx[8] = 1.0*reaction_16_k1;
    dwdx[9] = 1.0*reaction_32_k1;
    dwdx[10] = -1.0*reaction_15_V*species_11/pow(reaction_15_Km + species_11, 2) + 1.0*reaction_15_V/(reaction_15_Km + species_11);
    dwdx[11] = 1.0*reaction_31_k1;
    dwdx[12] = -1.0*reaction_47_V*species_11/pow(reaction_47_Km + species_11, 2) + 1.0*reaction_47_V/(reaction_47_Km + species_11);
    dwdx[13] = -1.0*reaction_16_k2;
    dwdx[14] = -1.0*reaction_21_V*species_12/pow(reaction_21_Km + species_12, 2) + 1.0*reaction_21_V/(reaction_21_Km + species_12);
    dwdx[15] = 1.0*reaction_23_k1*species_3;
    dwdx[16] = 1.0*reaction_46_k1;
    dwdx[17] = 1.0*reaction_19_k1;
    dwdx[18] = -1.0*reaction_37_V*species_13/pow(reaction_37_Km + species_13, 2) + 1.0*reaction_37_V/(reaction_37_Km + species_13);
    dwdx[19] = -1.0*reaction_42_V*species_13/pow(reaction_42_Km + species_13, 2) + 1.0*reaction_42_V/(reaction_42_Km + species_13);
    dwdx[20] = 1.0*reaction_17_k1;
    dwdx[21] = -1.0*reaction_33_V*species_14/pow(reaction_33_Km + species_14, 2) + 1.0*reaction_33_V/(reaction_33_Km + species_14);
    dwdx[22] = -1.0*reaction_45_V*species_14/pow(reaction_45_Km + species_14, 2) + 1.0*reaction_45_V/(reaction_45_Km + species_14);
    dwdx[23] = -1.0*reaction_23_k2;
    dwdx[24] = 1.0*reaction_38_k1;
    dwdx[25] = -1.0*reaction_39_V*species_15/pow(reaction_39_Km + species_15, 2) + 1.0*reaction_39_V/(reaction_39_Km + species_15);
    dwdx[26] = 1.0*reaction_18_k1;
    dwdx[27] = -1.0*reaction_35_V*species_2/pow(reaction_35_Km + species_2, 2) + 1.0*reaction_35_V/(reaction_35_Km + species_2);
    dwdx[28] = -1.0*reaction_41_V*species_2/pow(reaction_41_Km + species_2, 2) + 1.0*reaction_41_V/(reaction_41_Km + species_2);
    dwdx[29] = -1.0*pow(reaction_0_K, reaction_0_m)*reaction_0_m*reaction_0_vsb*pow(species_3, reaction_0_m)/(species_3*pow(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m), 2));
    dwdx[30] = -1.0*reaction_20_Vs*reaction_20_n*pow(species_3, 2*reaction_20_n)/(species_3*pow(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n), 2)) + 1.0*reaction_20_Vs*reaction_20_n*pow(species_3, reaction_20_n)/(species_3*(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n)));
    dwdx[31] = 1.0*reaction_23_k1*species_12;
    dwdx[32] = -1.0*reaction_36_V*species_3/pow(reaction_36_Km + species_3, 2) + 1.0*reaction_36_V/(reaction_36_Km + species_3);
    dwdx[33] = -1.0*reaction_4_k2;
    dwdx[34] = 1.0*reaction_40_k1;
    dwdx[35] = -1.0*reaction_9_Vs*reaction_9_n*pow(species_3, 2*reaction_9_n)/(species_3*pow(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n), 2)) + 1.0*reaction_9_Vs*reaction_9_n*pow(species_3, reaction_9_n)/(species_3*(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n)));
    dwdx[36] = 1.0*reaction_13_k1*species_8;
    dwdx[37] = 1.0*reaction_28_k1;
    dwdx[38] = -1.0*reaction_7_V*species_4/pow(reaction_7_Km + species_4, 2) + 1.0*reaction_7_V/(reaction_7_Km + species_4);
    dwdx[39] = -1.0*reaction_25_V*species_5/pow(reaction_25_Km + species_5, 2) + 1.0*reaction_25_V/(reaction_25_Km + species_5);
    dwdx[40] = 1.0*reaction_5_k;
    dwdx[41] = 1.0*reaction_6_k1;
    dwdx[42] = 1.0*reaction_30_k1;
    dwdx[43] = -1.0*reaction_43_V*species_6/pow(reaction_43_Km + species_6, 2) + 1.0*reaction_43_V/(reaction_43_Km + species_6);
    dwdx[44] = -1.0*reaction_8_V*species_6/pow(reaction_8_Km + species_6, 2) + 1.0*reaction_8_V/(reaction_8_Km + species_6);
    dwdx[45] = 1.0*reaction_10_k;
    dwdx[46] = 1.0*reaction_22_k1;
    dwdx[47] = -1.0*reaction_26_V*species_7/pow(reaction_26_Km + species_7, 2) + 1.0*reaction_26_V/(reaction_26_Km + species_7);
    dwdx[48] = -1.0*reaction_12_V*species_8/pow(reaction_12_Km + species_8, 2) + 1.0*reaction_12_V/(reaction_12_Km + species_8);
    dwdx[49] = 1.0*reaction_13_k1*species_4;
    dwdx[50] = 1.0*reaction_27_k1;
    dwdx[51] = -1.0*reaction_11_V*species_9/pow(reaction_11_Km + species_9, 2) + 1.0*reaction_11_V/(reaction_11_Km + species_9);
    dwdx[52] = 1.0*reaction_29_k1;
    dwdx[53] = -1.0*reaction_44_V*species_9/pow(reaction_44_Km + species_9, 2) + 1.0*reaction_44_V/(reaction_44_Km + species_9);
}