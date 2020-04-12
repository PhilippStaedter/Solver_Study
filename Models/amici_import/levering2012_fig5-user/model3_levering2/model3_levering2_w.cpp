#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model3_levering2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = reaction_1_Vmax*species_17*species_20*species_5/(reaction_1_KmS1*reaction_1_KmS2*reaction_1_KmS3*((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1));
    w[1] = reaction_10_KmI*reaction_10_Vmax*species_7*species_8*species_9/(reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1));
    w[2] = reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_7/(reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1));
    w[3] = 0.034799999999999998*reaction_12_V*pow(species_10, reaction_12_h)/(pow(reaction_12_Shalve, reaction_12_h) + pow(species_10, reaction_12_h));
    w[4] = reaction_13_KmI*reaction_13_Vmax*species_22*(species_22 - species_8)/((reaction_13_KmA + species_22)*(reaction_13_KmI + species_10)*(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP));
    w[5] = 0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15/(reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1));
    w[6] = reaction_15_KmI*reaction_15_Vmax*species_14/(reaction_15_KmS*(reaction_15_KmI + species_6)*(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP));
    w[7] = reaction_16_k*(-species_15 + species_20);
    w[8] = 0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_17*species_2/(reaction_17_KmS1*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1));
    w[9] = 0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*species_16*species_8/(reaction_18_KmS*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS));
    w[10] = 0.034799999999999998*reaction_19_Vmax*species_18*species_3/(reaction_19_KmS1*reaction_19_KmS2*((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1));
    w[11] = 0.034799999999999998*reaction_2_Vmax*species_1*species_10/(reaction_2_KmS1*reaction_2_KmS2*((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1));
    w[12] = 0.034799999999999998*reaction_20_Vmax*species_19/(reaction_20_KmS*(1 + species_19/reaction_20_KmS + species_18/reaction_20_KmP));
    w[13] = reaction_21_Vmax*species_10*species_22/(reaction_21_KmS1*reaction_21_KmS2*((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1));
    w[14] = 0.034799999999999998*(reaction_3_Vmax*species_2/reaction_3_KmS1 - reaction_3_Vmax*pow(species_3, 2)/(reaction_3_Keq*reaction_3_KmS1))/(1 + species_2/reaction_3_KmS1 + species_3/reaction_3_KmP1 + pow(species_3, 2)/pow(reaction_3_KmP1, 2));
    w[15] = 0.034799999999999998*reaction_4_Vmax*species_1*species_16/(reaction_4_KmS*(reaction_4_KmA + species_16)*((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS));
    w[16] = 0.034799999999999998*reaction_5_KmI*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/((reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1));
    w[17] = 0.034799999999999998*(reaction_6_Vmax*species_4*species_9/(reaction_6_KmS1*reaction_6_KmS2) - reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*reaction_6_KmS2))/((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1);
    w[18] = 0.034799999999999998*reaction_7_KmI*species_1*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/((reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1));
    w[19] = 0.034799999999999998*reaction_8_KmI*species_2*species_8*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/((reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1));
    w[20] = reaction_9_KmI*(reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/((reaction_9_KmI + species_3)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1));
}