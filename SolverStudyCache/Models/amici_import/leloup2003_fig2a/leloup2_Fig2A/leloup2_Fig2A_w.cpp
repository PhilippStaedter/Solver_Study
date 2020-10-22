#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_leloup2_Fig2A(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*pow(reaction_0_K, reaction_0_m)*reaction_0_vsb/(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m));
    w[1] = 1.0*reaction_1_k*species_0;
    w[2] = 1.0*reaction_10_k*species_7;
    w[3] = 1.0*reaction_11_V*species_9/(reaction_11_Km + species_9);
    w[4] = 1.0*reaction_12_V*species_8/(reaction_12_Km + species_8);
    w[5] = 1.0*reaction_13_k1*species_4*species_8 - 1.0*reaction_13_k2*species_10;
    w[6] = 1.0*reaction_14_V*species_10/(reaction_14_Km + species_10);
    w[7] = 1.0*reaction_15_V*species_11/(reaction_15_Km + species_11);
    w[8] = 1.0*reaction_16_k1*species_10 - 1.0*reaction_16_k2*species_12;
    w[9] = 1.0*reaction_17_k1*species_14;
    w[10] = 1.0*reaction_18_k1*species_2;
    w[11] = 1.0*reaction_19_k1*species_13;
    w[12] = 1.0*reaction_2_k1*species_0;
    w[13] = 1.0*reaction_20_Vs*pow(species_3, reaction_20_n)/(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n));
    w[14] = 1.0*reaction_21_V*species_12/(reaction_21_Km + species_12);
    w[15] = 1.0*reaction_22_k1*species_7;
    w[16] = 1.0*reaction_23_k1*species_12*species_3 - 1.0*reaction_23_k2*species_15;
    w[17] = 1.0*reaction_24_V*species_0/(reaction_24_Km + species_0);
    w[18] = 1.0*reaction_25_V*species_5/(reaction_25_Km + species_5);
    w[19] = 1.0*reaction_26_V*species_7/(reaction_26_Km + species_7);
    w[20] = 1.0*reaction_27_k1*species_8;
    w[21] = 1.0*reaction_28_k1*species_4;
    w[22] = 1.0*reaction_29_k1*species_9;
    w[23] = 1.0*reaction_3_V*species_1/(reaction_3_Km + species_1);
    w[24] = 1.0*reaction_30_k1*species_6;
    w[25] = 1.0*reaction_31_k1*species_11;
    w[26] = 1.0*reaction_32_k1*species_10;
    w[27] = 1.0*reaction_33_V*species_14/(reaction_33_Km + species_14);
    w[28] = 1.0*reaction_34_k1*species_1;
    w[29] = 1.0*reaction_35_V*species_2/(reaction_35_Km + species_2);
    w[30] = 1.0*reaction_36_V*species_3/(reaction_36_Km + species_3);
    w[31] = 1.0*reaction_37_V*species_13/(reaction_37_Km + species_13);
    w[32] = 1.0*reaction_38_k1*species_15;
    w[33] = 1.0*reaction_39_V*species_15/(reaction_39_Km + species_15);
    w[34] = 1.0*reaction_4_k1*species_1 - 1.0*reaction_4_k2*species_3;
    w[35] = 1.0*reaction_40_k1*species_3;
    w[36] = 1.0*reaction_41_V*species_2/(reaction_41_Km + species_2);
    w[37] = 1.0*reaction_42_V*species_13/(reaction_42_Km + species_13);
    w[38] = 1.0*reaction_43_V*species_6/(reaction_43_Km + species_6);
    w[39] = 1.0*reaction_44_V*species_9/(reaction_44_Km + species_9);
    w[40] = 1.0*reaction_45_V*species_14/(reaction_45_Km + species_14);
    w[41] = 1.0*reaction_46_k1*species_12;
    w[42] = 1.0*reaction_47_V*species_11/(reaction_47_Km + species_11);
    w[43] = 1.0*reaction_5_k*species_5;
    w[44] = 1.0*reaction_6_k1*species_5;
    w[45] = 1.0*reaction_7_V*species_4/(reaction_7_Km + species_4);
    w[46] = 1.0*reaction_8_V*species_6/(reaction_8_Km + species_6);
    w[47] = 1.0*reaction_9_Vs*pow(species_3, reaction_9_n)/(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n));
}