#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model3_levering2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 2:
            dwdp[0] = reaction_1_Vmax*species_1*species_17*species_20*species_5*(1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)/(pow(reaction_1_KmP3, 2)*reaction_1_KmS1*reaction_1_KmS2*reaction_1_KmS3*pow((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1, 2));
            break;
        case 3:
            dwdp[0] = reaction_1_Vmax*species_17*species_20*species_5*species_6*(1 + species_17/reaction_1_KmP1)*(1 + species_1/reaction_1_KmP3)/(pow(reaction_1_KmP2, 2)*reaction_1_KmS1*reaction_1_KmS2*reaction_1_KmS3*pow((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1, 2));
            break;
        case 4:
            dwdp[0] = reaction_1_Vmax*pow(species_17, 2)*species_20*species_5*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3)/(pow(reaction_1_KmP1, 2)*reaction_1_KmS1*reaction_1_KmS2*reaction_1_KmS3*pow((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1, 2));
            break;
        case 5:
            dwdp[0] = -reaction_1_Vmax*species_17*species_20*species_5/(reaction_1_KmS1*reaction_1_KmS2*pow(reaction_1_KmS3, 2)*((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1)) + reaction_1_Vmax*species_17*pow(species_20, 2)*species_5*(1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)/(reaction_1_KmS1*reaction_1_KmS2*pow(reaction_1_KmS3, 3)*pow((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1, 2));
            break;
        case 6:
            dwdp[0] = -reaction_1_Vmax*species_17*species_20*species_5/(reaction_1_KmS1*pow(reaction_1_KmS2, 2)*reaction_1_KmS3*((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1)) + reaction_1_Vmax*species_17*species_20*pow(species_5, 2)*(1 + species_17/reaction_1_KmS1)*(1 + species_20/reaction_1_KmS3)/(reaction_1_KmS1*pow(reaction_1_KmS2, 3)*reaction_1_KmS3*pow((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1, 2));
            break;
        case 7:
            dwdp[0] = -reaction_1_Vmax*species_17*species_20*species_5/(pow(reaction_1_KmS1, 2)*reaction_1_KmS2*reaction_1_KmS3*((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1)) + reaction_1_Vmax*pow(species_17, 2)*species_20*species_5*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3)/(pow(reaction_1_KmS1, 3)*reaction_1_KmS2*reaction_1_KmS3*pow((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1, 2));
            break;
        case 8:
            dwdp[0] = species_17*species_20*species_5/(reaction_1_KmS1*reaction_1_KmS2*reaction_1_KmS3*((1 + species_17/reaction_1_KmP1)*(1 + species_6/reaction_1_KmP2)*(1 + species_1/reaction_1_KmP3) + (1 + species_17/reaction_1_KmS1)*(1 + species_5/reaction_1_KmS2)*(1 + species_20/reaction_1_KmS3) - 1));
            break;
        case 9:
            dwdp[1] = reaction_10_KmI*reaction_10_Vmax*species_10*species_7*species_8*species_9*(1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)/(reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*pow(reaction_10_kmP3, 2)*(reaction_10_KmI + species_2)*pow((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1, 2));
            break;
        case 10:
            dwdp[1] = reaction_10_KmI*reaction_10_Vmax*species_13*species_7*species_8*species_9*(1 + species_25/reaction_10_KmP1)*(1 + species_10/reaction_10_kmP3)/(pow(reaction_10_KmP2, 2)*reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*pow((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1, 2));
            break;
        case 11:
            dwdp[1] = reaction_10_KmI*reaction_10_Vmax*species_25*species_7*species_8*species_9*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3)/(pow(reaction_10_KmP1, 2)*reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*pow((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1, 2));
            break;
        case 12:
            dwdp[1] = -reaction_10_KmI*reaction_10_Vmax*species_7*species_8*species_9/(reaction_10_KmS1*reaction_10_KmS2*pow(reaction_10_KmS3, 2)*(reaction_10_KmI + species_2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1)) + reaction_10_KmI*reaction_10_Vmax*species_7*species_8*pow(species_9, 2)*(1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)/(reaction_10_KmS1*reaction_10_KmS2*pow(reaction_10_KmS3, 3)*(reaction_10_KmI + species_2)*pow((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1, 2));
            break;
        case 13:
            dwdp[1] = -reaction_10_KmI*reaction_10_Vmax*species_7*species_8*species_9/(reaction_10_KmS1*pow(reaction_10_KmS2, 2)*reaction_10_KmS3*(reaction_10_KmI + species_2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1)) + reaction_10_KmI*reaction_10_Vmax*species_7*pow(species_8, 2)*species_9*(1 + species_7/reaction_10_KmS1)*(1 + species_9/reaction_10_KmS3)/(reaction_10_KmS1*pow(reaction_10_KmS2, 3)*reaction_10_KmS3*(reaction_10_KmI + species_2)*pow((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1, 2));
            break;
        case 14:
            dwdp[1] = -reaction_10_KmI*reaction_10_Vmax*species_7*species_8*species_9/(pow(reaction_10_KmS1, 2)*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1)) + reaction_10_KmI*reaction_10_Vmax*pow(species_7, 2)*species_8*species_9*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3)/(pow(reaction_10_KmS1, 3)*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*pow((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1, 2));
            break;
        case 15:
            dwdp[1] = reaction_10_KmI*species_7*species_8*species_9/(reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1));
            break;
        case 16:
            dwdp[1] = -reaction_10_KmI*reaction_10_Vmax*species_7*species_8*species_9/(reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*pow(reaction_10_KmI + species_2, 2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1)) + reaction_10_Vmax*species_7*species_8*species_9/(reaction_10_KmS1*reaction_10_KmS2*reaction_10_KmS3*(reaction_10_KmI + species_2)*((1 + species_25/reaction_10_KmP1)*(1 + species_13/reaction_10_KmP2)*(1 + species_10/reaction_10_kmP3) + (1 + species_7/reaction_10_KmS1)*(1 + species_8/reaction_10_KmS2)*(1 + species_9/reaction_10_KmS3) - 1));
            break;
        case 17:
            dwdp[2] = reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_13*species_7*(1 + species_24/reaction_11_KmP1)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2))/(pow(reaction_11_KmP4, 2)*reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*pow((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1, 2));
            break;
        case 18:
            dwdp[2] = -reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_7*(1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(-species_11/pow(reaction_11_KmP2, 2) - 2*pow(species_11, 2)/pow(reaction_11_KmP2, 3))/(reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*pow((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1, 2));
            break;
        case 19:
            dwdp[2] = reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_24*species_7*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2))/(pow(reaction_11_KmP1, 2)*reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*pow((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1, 2));
            break;
        case 20:
            dwdp[2] = -reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_7*(1 + species_7/reaction_11_KmS1)*(-species_12/pow(reaction_11_KmS2, 2) - 2*pow(species_12, 2)/pow(reaction_11_KmS2, 3))/(reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*pow((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1, 2)) - 2*reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_7/(reaction_11_KmS1*pow(reaction_11_KmS2, 3)*(reaction_11_KmI + species_10)*((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1));
            break;
        case 21:
            dwdp[2] = -reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_7/(pow(reaction_11_KmS1, 2)*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1)) + reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*pow(species_7, 2)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2))/(pow(reaction_11_KmS1, 3)*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*pow((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1, 2));
            break;
        case 22:
            dwdp[2] = reaction_11_KmI*pow(species_12, 2)*species_7/(reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1));
            break;
        case 23:
            dwdp[2] = -reaction_11_KmI*reaction_11_Vmax*pow(species_12, 2)*species_7/(reaction_11_KmS1*pow(reaction_11_KmS2, 2)*pow(reaction_11_KmI + species_10, 2)*((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1)) + reaction_11_Vmax*pow(species_12, 2)*species_7/(reaction_11_KmS1*pow(reaction_11_KmS2, 2)*(reaction_11_KmI + species_10)*((1 + species_24/reaction_11_KmP1)*(1 + species_13/reaction_11_KmP4)*(1 + species_11/reaction_11_KmP2 + pow(species_11, 2)/pow(reaction_11_KmP2, 2)) + (1 + species_7/reaction_11_KmS1)*(1 + species_12/reaction_11_KmS2 + pow(species_12, 2)/pow(reaction_11_KmS2, 2)) - 1));
            break;
        case 24:
            dwdp[3] = 0.034799999999999998*reaction_12_V*pow(species_10, reaction_12_h)*log(species_10)/(pow(reaction_12_Shalve, reaction_12_h) + pow(species_10, reaction_12_h)) + 0.034799999999999998*reaction_12_V*pow(species_10, reaction_12_h)*(-pow(reaction_12_Shalve, reaction_12_h)*log(reaction_12_Shalve) - pow(species_10, reaction_12_h)*log(species_10))/pow(pow(reaction_12_Shalve, reaction_12_h) + pow(species_10, reaction_12_h), 2);
            break;
        case 25:
            dwdp[3] = 0.034799999999999998*pow(species_10, reaction_12_h)/(pow(reaction_12_Shalve, reaction_12_h) + pow(species_10, reaction_12_h));
            break;
        case 26:
            dwdp[3] = -0.034799999999999998*pow(reaction_12_Shalve, reaction_12_h)*reaction_12_V*reaction_12_h*pow(species_10, reaction_12_h)/(reaction_12_Shalve*pow(pow(reaction_12_Shalve, reaction_12_h) + pow(species_10, reaction_12_h), 2));
            break;
        case 27:
            dwdp[4] = reaction_13_KmI*reaction_13_Vmax*species_22*species_8*(species_22 - species_8)/(pow(reaction_13_KmP, 2)*(reaction_13_KmA + species_22)*(reaction_13_KmI + species_10)*pow(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP, 2));
            break;
        case 28:
            dwdp[4] = reaction_13_KmI*reaction_13_Vmax*pow(species_22, 2)*(species_22 - species_8)/(pow(reaction_13_KmS, 2)*(reaction_13_KmA + species_22)*(reaction_13_KmI + species_10)*pow(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP, 2));
            break;
        case 29:
            dwdp[4] = reaction_13_KmI*species_22*(species_22 - species_8)/((reaction_13_KmA + species_22)*(reaction_13_KmI + species_10)*(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP));
            break;
        case 30:
            dwdp[4] = -reaction_13_KmI*reaction_13_Vmax*species_22*(species_22 - species_8)/(pow(reaction_13_KmA + species_22, 2)*(reaction_13_KmI + species_10)*(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP));
            break;
        case 31:
            dwdp[4] = -reaction_13_KmI*reaction_13_Vmax*species_22*(species_22 - species_8)/((reaction_13_KmA + species_22)*pow(reaction_13_KmI + species_10, 2)*(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP)) + reaction_13_Vmax*species_22*(species_22 - species_8)/((reaction_13_KmA + species_22)*(reaction_13_KmI + species_10)*(1 + species_22/reaction_13_KmS + species_8/reaction_13_KmP));
            break;
        case 32:
            dwdp[5] = 0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15*species_9*(1 + species_1/reaction_14_KmP1)/(pow(reaction_14_KmP2, 2)*reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*pow((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1, 2));
            break;
        case 33:
            dwdp[5] = 0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_1*species_10*species_15*(1 + species_9/reaction_14_KmP2)/(pow(reaction_14_KmP1, 2)*reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*pow((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1, 2));
            break;
        case 34:
            dwdp[5] = -0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15/(reaction_14_KmS1*pow(reaction_14_KmS2, 2)*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1)) + 0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*pow(species_10, 2)*species_15*(1 + species_15/reaction_14_KmS1)/(reaction_14_KmS1*pow(reaction_14_KmS2, 3)*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*pow((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1, 2));
            break;
        case 35:
            dwdp[5] = -0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15/(pow(reaction_14_KmS1, 2)*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1)) + 0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*pow(species_15, 2)*(1 + species_10/reaction_14_KmS2)/(pow(reaction_14_KmS1, 3)*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*pow((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1, 2));
            break;
        case 36:
            dwdp[5] = 0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*species_10*species_15/(reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1));
            break;
        case 37:
            dwdp[5] = -0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15/(reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*pow(reaction_14_KmI2 + species_9, 2)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1)) + 0.034799999999999998*reaction_14_KmI1*reaction_14_Vmax*species_10*species_15/(reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1));
            break;
        case 38:
            dwdp[5] = -0.034799999999999998*reaction_14_KmI1*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15/(reaction_14_KmS1*reaction_14_KmS2*pow(reaction_14_KmI1 + species_1, 2)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1)) + 0.034799999999999998*reaction_14_KmI2*reaction_14_Vmax*species_10*species_15/(reaction_14_KmS1*reaction_14_KmS2*(reaction_14_KmI1 + species_1)*(reaction_14_KmI2 + species_9)*((1 + species_1/reaction_14_KmP1)*(1 + species_9/reaction_14_KmP2) + (1 + species_15/reaction_14_KmS1)*(1 + species_10/reaction_14_KmS2) - 1));
            break;
        case 39:
            dwdp[6] = reaction_15_KmI*reaction_15_Vmax*species_14*species_21/(pow(reaction_15_KmP, 2)*reaction_15_KmS*(reaction_15_KmI + species_6)*pow(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP, 2));
            break;
        case 40:
            dwdp[6] = -reaction_15_KmI*reaction_15_Vmax*species_14/(pow(reaction_15_KmS, 2)*(reaction_15_KmI + species_6)*(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP)) + reaction_15_KmI*reaction_15_Vmax*pow(species_14, 2)/(pow(reaction_15_KmS, 3)*(reaction_15_KmI + species_6)*pow(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP, 2));
            break;
        case 41:
            dwdp[6] = reaction_15_KmI*species_14/(reaction_15_KmS*(reaction_15_KmI + species_6)*(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP));
            break;
        case 42:
            dwdp[6] = -reaction_15_KmI*reaction_15_Vmax*species_14/(reaction_15_KmS*pow(reaction_15_KmI + species_6, 2)*(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP)) + reaction_15_Vmax*species_14/(reaction_15_KmS*(reaction_15_KmI + species_6)*(1 + species_14/reaction_15_KmS + species_21/reaction_15_KmP));
            break;
        case 43:
            dwdp[7] = -species_15 + species_20;
            break;
        case 44:
            dwdp[8] = 0.034799999999999998*reaction_17_KmI*species_10*species_17*species_2/(reaction_17_KmS1*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1));
            break;
        case 45:
            dwdp[8] = -0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_17*species_2/(reaction_17_KmS1*reaction_17_KmS2*pow(reaction_17_KmA + species_2, 2)*(reaction_17_KmI + species_8)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1));
            break;
        case 46:
            dwdp[8] = -0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_17*species_2/(reaction_17_KmS1*reaction_17_KmS2*(reaction_17_KmA + species_2)*pow(reaction_17_KmI + species_8, 2)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1)) + 0.034799999999999998*reaction_17_Vmax*species_10*species_17*species_2/(reaction_17_KmS1*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1));
            break;
        case 47:
            dwdp[8] = 0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_17*species_2*species_9*(1 + species_16/reaction_17_KmP1)/(pow(reaction_17_KmP2, 2)*reaction_17_KmS1*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*pow((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1, 2));
            break;
        case 48:
            dwdp[8] = 0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_16*species_17*species_2*(1 + species_9/reaction_17_KmP2)/(pow(reaction_17_KmP1, 2)*reaction_17_KmS1*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*pow((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1, 2));
            break;
        case 49:
            dwdp[8] = -0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_17*species_2/(reaction_17_KmS1*pow(reaction_17_KmS2, 2)*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1)) + 0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*pow(species_10, 2)*species_17*species_2*(1 + species_17/reaction_17_KmS1)/(reaction_17_KmS1*pow(reaction_17_KmS2, 3)*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*pow((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1, 2));
            break;
        case 50:
            dwdp[8] = -0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*species_17*species_2/(pow(reaction_17_KmS1, 2)*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1)) + 0.034799999999999998*reaction_17_KmI*reaction_17_Vmax*species_10*pow(species_17, 2)*species_2*(1 + species_10/reaction_17_KmS2)/(pow(reaction_17_KmS1, 3)*reaction_17_KmS2*(reaction_17_KmA + species_2)*(reaction_17_KmI + species_8)*pow((1 + species_16/reaction_17_KmP1)*(1 + species_9/reaction_17_KmP2) + (1 + species_17/reaction_17_KmS1)*(1 + species_10/reaction_17_KmS2) - 1, 2));
            break;
        case 51:
            dwdp[9] = 0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*species_16*pow(species_8, 2)*(1 + species_17/reaction_18_KmP1)/(pow(reaction_18_KmP2, 2)*reaction_18_KmS*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*pow((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS, 2));
            break;
        case 52:
            dwdp[9] = 0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*species_16*species_17*species_8*(1 + species_8/reaction_18_KmP2)/(pow(reaction_18_KmP1, 2)*reaction_18_KmS*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*pow((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS, 2));
            break;
        case 53:
            dwdp[9] = -0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*species_16*species_8/(pow(reaction_18_KmS, 2)*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS)) + 0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*pow(species_16, 2)*species_8/(pow(reaction_18_KmS, 3)*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*pow((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS, 2));
            break;
        case 54:
            dwdp[9] = 0.034799999999999998*reaction_18_KmI*species_16*species_8/(reaction_18_KmS*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS));
            break;
        case 55:
            dwdp[9] = -0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*species_16*species_8/(reaction_18_KmS*pow(reaction_18_KmA + species_8, 2)*(reaction_18_KmI + species_10)*((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS));
            break;
        case 56:
            dwdp[9] = -0.034799999999999998*reaction_18_KmI*reaction_18_Vmax*species_16*species_8/(reaction_18_KmS*(reaction_18_KmA + species_8)*pow(reaction_18_KmI + species_10, 2)*((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS)) + 0.034799999999999998*reaction_18_Vmax*species_16*species_8/(reaction_18_KmS*(reaction_18_KmA + species_8)*(reaction_18_KmI + species_10)*((1 + species_17/reaction_18_KmP1)*(1 + species_8/reaction_18_KmP2) + species_16/reaction_18_KmS));
            break;
        case 57:
            dwdp[10] = 0.034799999999999998*species_18*species_3/(reaction_19_KmS1*reaction_19_KmS2*((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1));
            break;
        case 58:
            dwdp[10] = 0.034799999999999998*reaction_19_Vmax*species_18*species_19*species_3*(1 + species_5/reaction_19_KmP1)/(pow(reaction_19_KmP2, 2)*reaction_19_KmS1*reaction_19_KmS2*pow((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1, 2));
            break;
        case 59:
            dwdp[10] = 0.034799999999999998*reaction_19_Vmax*species_18*species_3*species_5*(1 + species_19/reaction_19_KmP2)/(pow(reaction_19_KmP1, 2)*reaction_19_KmS1*reaction_19_KmS2*pow((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1, 2));
            break;
        case 60:
            dwdp[10] = -0.034799999999999998*reaction_19_Vmax*species_18*species_3/(reaction_19_KmS1*pow(reaction_19_KmS2, 2)*((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1)) + 0.034799999999999998*reaction_19_Vmax*pow(species_18, 2)*species_3*(1 + species_3/reaction_19_KmS1)/(reaction_19_KmS1*pow(reaction_19_KmS2, 3)*pow((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1, 2));
            break;
        case 61:
            dwdp[10] = -0.034799999999999998*reaction_19_Vmax*species_18*species_3/(pow(reaction_19_KmS1, 2)*reaction_19_KmS2*((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1)) + 0.034799999999999998*reaction_19_Vmax*species_18*pow(species_3, 2)*(1 + species_18/reaction_19_KmS2)/(pow(reaction_19_KmS1, 3)*reaction_19_KmS2*pow((1 + species_5/reaction_19_KmP1)*(1 + species_19/reaction_19_KmP2) + (1 + species_3/reaction_19_KmS1)*(1 + species_18/reaction_19_KmS2) - 1, 2));
            break;
        case 62:
            dwdp[11] = 0.034799999999999998*species_1*species_10/(reaction_2_KmS1*reaction_2_KmS2*((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1));
            break;
        case 63:
            dwdp[11] = 0.034799999999999998*reaction_2_Vmax*species_1*species_10*species_9*(1 + species_2/reaction_2_KmP1)/(pow(reaction_2_KmP2, 2)*reaction_2_KmS1*reaction_2_KmS2*pow((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1, 2));
            break;
        case 64:
            dwdp[11] = 0.034799999999999998*reaction_2_Vmax*species_1*species_10*species_2*(1 + species_9/reaction_2_KmP2)/(pow(reaction_2_KmP1, 2)*reaction_2_KmS1*reaction_2_KmS2*pow((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1, 2));
            break;
        case 65:
            dwdp[11] = -0.034799999999999998*reaction_2_Vmax*species_1*species_10/(reaction_2_KmS1*pow(reaction_2_KmS2, 2)*((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1)) + 0.034799999999999998*reaction_2_Vmax*species_1*pow(species_10, 2)*(1 + species_1/reaction_2_KmS1)/(reaction_2_KmS1*pow(reaction_2_KmS2, 3)*pow((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1, 2));
            break;
        case 66:
            dwdp[11] = -0.034799999999999998*reaction_2_Vmax*species_1*species_10/(pow(reaction_2_KmS1, 2)*reaction_2_KmS2*((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1)) + 0.034799999999999998*reaction_2_Vmax*pow(species_1, 2)*species_10*(1 + species_10/reaction_2_KmS2)/(pow(reaction_2_KmS1, 3)*reaction_2_KmS2*pow((1 + species_2/reaction_2_KmP1)*(1 + species_9/reaction_2_KmP2) + (1 + species_1/reaction_2_KmS1)*(1 + species_10/reaction_2_KmS2) - 1, 2));
            break;
        case 67:
            dwdp[12] = 0.034799999999999998*reaction_20_Vmax*species_18*species_19/(pow(reaction_20_KmP, 2)*reaction_20_KmS*pow(1 + species_19/reaction_20_KmS + species_18/reaction_20_KmP, 2));
            break;
        case 68:
            dwdp[12] = -0.034799999999999998*reaction_20_Vmax*species_19/(pow(reaction_20_KmS, 2)*(1 + species_19/reaction_20_KmS + species_18/reaction_20_KmP)) + 0.034799999999999998*reaction_20_Vmax*pow(species_19, 2)/(pow(reaction_20_KmS, 3)*pow(1 + species_19/reaction_20_KmS + species_18/reaction_20_KmP, 2));
            break;
        case 69:
            dwdp[12] = 0.034799999999999998*species_19/(reaction_20_KmS*(1 + species_19/reaction_20_KmS + species_18/reaction_20_KmP));
            break;
        case 70:
            dwdp[13] = -reaction_21_Vmax*species_10*species_22*(1 + species_9/reaction_21_KmP1)*(-species_8/pow(reaction_21_KmP2, 2) - 2*pow(species_8, 2)/pow(reaction_21_KmP2, 3))/(reaction_21_KmS1*reaction_21_KmS2*pow((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1, 2));
            break;
        case 71:
            dwdp[13] = reaction_21_Vmax*species_10*species_22*species_9*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2))/(pow(reaction_21_KmP1, 2)*reaction_21_KmS1*reaction_21_KmS2*pow((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1, 2));
            break;
        case 72:
            dwdp[13] = -reaction_21_Vmax*species_10*species_22/(reaction_21_KmS1*pow(reaction_21_KmS2, 2)*((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1)) + reaction_21_Vmax*pow(species_10, 2)*species_22*(1 + species_22/reaction_21_KmS1)/(reaction_21_KmS1*pow(reaction_21_KmS2, 3)*pow((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1, 2));
            break;
        case 73:
            dwdp[13] = -reaction_21_Vmax*species_10*species_22/(pow(reaction_21_KmS1, 2)*reaction_21_KmS2*((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1)) + reaction_21_Vmax*species_10*pow(species_22, 2)*(1 + species_10/reaction_21_KmS2)/(pow(reaction_21_KmS1, 3)*reaction_21_KmS2*pow((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1, 2));
            break;
        case 74:
            dwdp[13] = species_10*species_22/(reaction_21_KmS1*reaction_21_KmS2*((1 + species_9/reaction_21_KmP1)*(1 + species_8/reaction_21_KmP2 + pow(species_8, 2)/pow(reaction_21_KmP2, 2)) + (1 + species_22/reaction_21_KmS1)*(1 + species_10/reaction_21_KmS2) - 1));
            break;
        case 75:
            dwdp[14] = (species_3/pow(reaction_3_KmP1, 2) + 2*pow(species_3, 2)/pow(reaction_3_KmP1, 3))*(0.034799999999999998*reaction_3_Vmax*species_2/reaction_3_KmS1 - 0.034799999999999998*reaction_3_Vmax*pow(species_3, 2)/(reaction_3_Keq*reaction_3_KmS1))/pow(1 + species_2/reaction_3_KmS1 + species_3/reaction_3_KmP1 + pow(species_3, 2)/pow(reaction_3_KmP1, 2), 2);
            break;
        case 76:
            dwdp[14] = 0.034799999999999998*reaction_3_Vmax*pow(species_3, 2)/(pow(reaction_3_Keq, 2)*reaction_3_KmS1*(1 + species_2/reaction_3_KmS1 + species_3/reaction_3_KmP1 + pow(species_3, 2)/pow(reaction_3_KmP1, 2)));
            break;
        case 77:
            dwdp[14] = 0.034799999999999998*(-reaction_3_Vmax*species_2/pow(reaction_3_KmS1, 2) + reaction_3_Vmax*pow(species_3, 2)/(reaction_3_Keq*pow(reaction_3_KmS1, 2)))/(1 + species_2/reaction_3_KmS1 + species_3/reaction_3_KmP1 + pow(species_3, 2)/pow(reaction_3_KmP1, 2)) + 0.034799999999999998*species_2*(reaction_3_Vmax*species_2/reaction_3_KmS1 - reaction_3_Vmax*pow(species_3, 2)/(reaction_3_Keq*reaction_3_KmS1))/(pow(reaction_3_KmS1, 2)*pow(1 + species_2/reaction_3_KmS1 + species_3/reaction_3_KmP1 + pow(species_3, 2)/pow(reaction_3_KmP1, 2), 2));
            break;
        case 78:
            dwdp[14] = 0.034799999999999998*(species_2/reaction_3_KmS1 - pow(species_3, 2)/(reaction_3_Keq*reaction_3_KmS1))/(1 + species_2/reaction_3_KmS1 + species_3/reaction_3_KmP1 + pow(species_3, 2)/pow(reaction_3_KmP1, 2));
            break;
        case 79:
            dwdp[15] = 0.034799999999999998*species_1*species_16/(reaction_4_KmS*(reaction_4_KmA + species_16)*((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS));
            break;
        case 80:
            dwdp[15] = -0.034799999999999998*reaction_4_Vmax*species_1*species_16/(reaction_4_KmS*pow(reaction_4_KmA + species_16, 2)*((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS));
            break;
        case 81:
            dwdp[15] = 0.034799999999999998*reaction_4_Vmax*species_1*species_16*species_8*(1 + species_15/reaction_4_KmP1)/(pow(reaction_4_KmP2, 2)*reaction_4_KmS*(reaction_4_KmA + species_16)*pow((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS, 2));
            break;
        case 82:
            dwdp[15] = 0.034799999999999998*reaction_4_Vmax*species_1*species_15*species_16*(1 + species_8/reaction_4_KmP2)/(pow(reaction_4_KmP1, 2)*reaction_4_KmS*(reaction_4_KmA + species_16)*pow((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS, 2));
            break;
        case 83:
            dwdp[15] = -0.034799999999999998*reaction_4_Vmax*species_1*species_16/(pow(reaction_4_KmS, 2)*(reaction_4_KmA + species_16)*((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS)) + 0.034799999999999998*reaction_4_Vmax*pow(species_1, 2)*species_16/(pow(reaction_4_KmS, 3)*(reaction_4_KmA + species_16)*pow((1 + species_15/reaction_4_KmP1)*(1 + species_8/reaction_4_KmP2) + species_1/reaction_4_KmS, 2));
            break;
        case 84:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*species_12*(1 + species_4/reaction_5_KmP1)*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/(pow(reaction_5_KmP2, 2)*(reaction_5_KmI + species_12)*pow((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1, 2));
            break;
        case 85:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*species_4*(1 + species_12/reaction_5_KmP2)*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/(pow(reaction_5_KmP1, 2)*(reaction_5_KmI + species_12)*pow((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1, 2));
            break;
        case 86:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*reaction_5_Vmax*species_12*species_4/(pow(reaction_5_Keq, 2)*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3*(reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1));
            break;
        case 87:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*(-reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*pow(reaction_5_KmS3, 2)) + reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*pow(reaction_5_KmS3, 2)))/((reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1)) + 0.034799999999999998*reaction_5_KmI*species_8*(1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/(pow(reaction_5_KmS3, 2)*(reaction_5_KmI + species_12)*pow((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1, 2));
            break;
        case 88:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*(-reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*pow(reaction_5_KmS2, 2)*reaction_5_KmS3) + reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*pow(reaction_5_KmS2, 2)*reaction_5_KmS3))/((reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1)) + 0.034799999999999998*reaction_5_KmI*species_11*(1 + species_3/reaction_5_KmS1)*(1 + species_8/reaction_5_KmS3)*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/(pow(reaction_5_KmS2, 2)*(reaction_5_KmI + species_12)*pow((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1, 2));
            break;
        case 89:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*(-reaction_5_Vmax*species_11*species_3*species_8/(pow(reaction_5_KmS1, 2)*reaction_5_KmS2*reaction_5_KmS3) + reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*pow(reaction_5_KmS1, 2)*reaction_5_KmS2*reaction_5_KmS3))/((reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1)) + 0.034799999999999998*reaction_5_KmI*species_3*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3)*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/(pow(reaction_5_KmS1, 2)*(reaction_5_KmI + species_12)*pow((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1, 2));
            break;
        case 90:
            dwdp[16] = 0.034799999999999998*reaction_5_KmI*(species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/((reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1));
            break;
        case 91:
            dwdp[16] = -0.034799999999999998*reaction_5_KmI*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/(pow(reaction_5_KmI + species_12, 2)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1)) + 0.034799999999999998*(reaction_5_Vmax*species_11*species_3*species_8/(reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3) - reaction_5_Vmax*species_12*species_4/(reaction_5_Keq*reaction_5_KmS1*reaction_5_KmS2*reaction_5_KmS3))/((reaction_5_KmI + species_12)*((1 + species_4/reaction_5_KmP1)*(1 + species_12/reaction_5_KmP2) + (1 + species_3/reaction_5_KmS1)*(1 + species_11/reaction_5_KmS2)*(1 + species_8/reaction_5_KmS3) - 1));
            break;
        case 92:
            dwdp[17] = species_10*(1 + species_5/reaction_6_KmP1)*(0.034799999999999998*reaction_6_Vmax*species_4*species_9/(reaction_6_KmS1*reaction_6_KmS2) - 0.034799999999999998*reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*reaction_6_KmS2))/(pow(reaction_6_KmP2, 2)*pow((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1, 2));
            break;
        case 93:
            dwdp[17] = species_5*(1 + species_10/reaction_6_KmP2)*(0.034799999999999998*reaction_6_Vmax*species_4*species_9/(reaction_6_KmS1*reaction_6_KmS2) - 0.034799999999999998*reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*reaction_6_KmS2))/(pow(reaction_6_KmP1, 2)*pow((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1, 2));
            break;
        case 94:
            dwdp[17] = 0.034799999999999998*reaction_6_Vmax*species_10*species_5/(pow(reaction_6_Keq, 2)*reaction_6_KmS1*reaction_6_KmS2*((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1));
            break;
        case 95:
            dwdp[17] = 0.034799999999999998*(-reaction_6_Vmax*species_4*species_9/(reaction_6_KmS1*pow(reaction_6_KmS2, 2)) + reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*pow(reaction_6_KmS2, 2)))/((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1) + 0.034799999999999998*species_9*(1 + species_4/reaction_6_KmS1)*(reaction_6_Vmax*species_4*species_9/(reaction_6_KmS1*reaction_6_KmS2) - reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*reaction_6_KmS2))/(pow(reaction_6_KmS2, 2)*pow((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1, 2));
            break;
        case 96:
            dwdp[17] = 0.034799999999999998*(-reaction_6_Vmax*species_4*species_9/(pow(reaction_6_KmS1, 2)*reaction_6_KmS2) + reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*pow(reaction_6_KmS1, 2)*reaction_6_KmS2))/((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1) + 0.034799999999999998*species_4*(1 + species_9/reaction_6_KmS2)*(reaction_6_Vmax*species_4*species_9/(reaction_6_KmS1*reaction_6_KmS2) - reaction_6_Vmax*species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*reaction_6_KmS2))/(pow(reaction_6_KmS1, 2)*pow((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1, 2));
            break;
        case 97:
            dwdp[17] = 0.034799999999999998*(species_4*species_9/(reaction_6_KmS1*reaction_6_KmS2) - species_10*species_5/(reaction_6_Keq*reaction_6_KmS1*reaction_6_KmS2))/((1 + species_5/reaction_6_KmP1)*(1 + species_10/reaction_6_KmP2) + (1 + species_4/reaction_6_KmS1)*(1 + species_9/reaction_6_KmS2) - 1);
            break;
        case 98:
            dwdp[18] = 0.034799999999999998*reaction_7_KmI*species_1*species_10*(1 + species_6/reaction_7_KmP1)*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/(pow(reaction_7_KmP2, 2)*(reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*pow((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1, 2));
            break;
        case 99:
            dwdp[18] = 0.034799999999999998*reaction_7_KmI*species_1*species_6*(1 + species_10/reaction_7_KmP2)*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/(pow(reaction_7_KmP1, 2)*(reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*pow((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1, 2));
            break;
        case 100:
            dwdp[18] = 0.034799999999999998*reaction_7_KmI*reaction_7_Vmax*species_1*species_10*species_6/(pow(reaction_7_Keq, 2)*reaction_7_KmS1*reaction_7_KmS2*(reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1));
            break;
        case 101:
            dwdp[18] = 0.034799999999999998*reaction_7_KmI*species_1*(-reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*pow(reaction_7_KmS2, 2)) + reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*pow(reaction_7_KmS2, 2)))/((reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1)) + 0.034799999999999998*reaction_7_KmI*species_1*species_9*(1 + species_5/reaction_7_KmS1)*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/(pow(reaction_7_KmS2, 2)*(reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*pow((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1, 2));
            break;
        case 102:
            dwdp[18] = 0.034799999999999998*reaction_7_KmI*species_1*(-reaction_7_Vmax*species_5*species_9/(pow(reaction_7_KmS1, 2)*reaction_7_KmS2) + reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*pow(reaction_7_KmS1, 2)*reaction_7_KmS2))/((reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1)) + 0.034799999999999998*reaction_7_KmI*species_1*species_5*(1 + species_9/reaction_7_KmS2)*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/(pow(reaction_7_KmS1, 2)*(reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*pow((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1, 2));
            break;
        case 103:
            dwdp[18] = 0.034799999999999998*reaction_7_KmI*species_1*(species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/((reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1));
            break;
        case 104:
            dwdp[18] = -0.034799999999999998*reaction_7_KmI*species_1*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/(pow(reaction_7_KmA + species_1, 2)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1));
            break;
        case 105:
            dwdp[18] = -0.034799999999999998*reaction_7_KmI*species_1*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/((reaction_7_KmA + species_1)*pow(reaction_7_KmI + species_8, 2)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1)) + 0.034799999999999998*species_1*(reaction_7_Vmax*species_5*species_9/(reaction_7_KmS1*reaction_7_KmS2) - reaction_7_Vmax*species_10*species_6/(reaction_7_Keq*reaction_7_KmS1*reaction_7_KmS2))/((reaction_7_KmA + species_1)*(reaction_7_KmI + species_8)*((1 + species_6/reaction_7_KmP1)*(1 + species_10/reaction_7_KmP2) + (1 + species_5/reaction_7_KmS1)*(1 + species_9/reaction_7_KmS2) - 1));
            break;
        case 106:
            dwdp[19] = 0.034799999999999998*reaction_8_KmI*species_11*species_2*species_8*(1 + species_14/reaction_8_KmP1)*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/(pow(reaction_8_KmP2, 2)*(reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*pow((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1, 2));
            break;
        case 107:
            dwdp[19] = 0.034799999999999998*reaction_8_KmI*species_14*species_2*species_8*(1 + species_11/reaction_8_KmP2)*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/(pow(reaction_8_KmP1, 2)*(reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*pow((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1, 2));
            break;
        case 108:
            dwdp[19] = 0.034799999999999998*reaction_8_KmI*reaction_8_Vmax*species_11*species_14*species_2*species_8/(pow(reaction_8_Keq, 2)*reaction_8_KmS1*reaction_8_KmS2*(reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1));
            break;
        case 109:
            dwdp[19] = 0.034799999999999998*reaction_8_KmI*species_2*species_8*(-reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*pow(reaction_8_KmS2, 2)) + reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*pow(reaction_8_KmS2, 2)))/((reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1)) + 0.034799999999999998*reaction_8_KmI*species_12*species_2*species_8*(1 + species_6/reaction_8_KmS1)*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/(pow(reaction_8_KmS2, 2)*(reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*pow((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1, 2));
            break;
        case 110:
            dwdp[19] = 0.034799999999999998*reaction_8_KmI*species_2*species_8*(-reaction_8_Vmax*species_12*species_6/(pow(reaction_8_KmS1, 2)*reaction_8_KmS2) + reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*pow(reaction_8_KmS1, 2)*reaction_8_KmS2))/((reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1)) + 0.034799999999999998*reaction_8_KmI*species_2*species_6*species_8*(1 + species_12/reaction_8_KmS2)*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/(pow(reaction_8_KmS1, 2)*(reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*pow((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1, 2));
            break;
        case 111:
            dwdp[19] = 0.034799999999999998*reaction_8_KmI*species_2*species_8*(species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/((reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1));
            break;
        case 112:
            dwdp[19] = -0.034799999999999998*reaction_8_KmI*species_2*species_8*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/((reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*pow(reaction_8_KmI + species_11, 2)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1)) + 0.034799999999999998*species_2*species_8*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/((reaction_8_KmA1 + species_2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1));
            break;
        case 113:
            dwdp[19] = -0.034799999999999998*reaction_8_KmI*species_2*species_8*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/((reaction_8_KmA1 + species_2)*pow(reaction_8_KmA2 + species_8, 2)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1));
            break;
        case 114:
            dwdp[19] = -0.034799999999999998*reaction_8_KmI*species_2*species_8*(reaction_8_Vmax*species_12*species_6/(reaction_8_KmS1*reaction_8_KmS2) - reaction_8_Vmax*species_11*species_14/(reaction_8_Keq*reaction_8_KmS1*reaction_8_KmS2))/(pow(reaction_8_KmA1 + species_2, 2)*(reaction_8_KmA2 + species_8)*(reaction_8_KmI + species_11)*((1 + species_14/reaction_8_KmP1)*(1 + species_11/reaction_8_KmP2) + (1 + species_6/reaction_8_KmS1)*(1 + species_12/reaction_8_KmS2) - 1));
            break;
        case 115:
            dwdp[20] = reaction_9_KmI*species_23*(1 + species_7/reaction_9_KmP1)*(reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/(pow(reaction_9_KmP2, 2)*(reaction_9_KmI + species_3)*pow((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1, 2));
            break;
        case 116:
            dwdp[20] = reaction_9_KmI*species_7*(1 + species_23/reaction_9_KmP2)*(reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/(pow(reaction_9_KmP1, 2)*(reaction_9_KmI + species_3)*pow((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1, 2));
            break;
        case 117:
            dwdp[20] = reaction_9_KmI*reaction_9_Vmax*species_23*species_7/(pow(reaction_9_Keq, 2)*reaction_9_KmS1*reaction_9_KmS2*(reaction_9_KmI + species_3)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1));
            break;
        case 118:
            dwdp[20] = reaction_9_KmI*(-reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*pow(reaction_9_KmS2, 2)) + reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*pow(reaction_9_KmS2, 2)))/((reaction_9_KmI + species_3)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1)) + reaction_9_KmI*species_13*(1 + species_6/reaction_9_KmS1)*(reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/(pow(reaction_9_KmS2, 2)*(reaction_9_KmI + species_3)*pow((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1, 2));
            break;
        case 119:
            dwdp[20] = reaction_9_KmI*(-reaction_9_Vmax*species_13*species_6/(pow(reaction_9_KmS1, 2)*reaction_9_KmS2) + reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*pow(reaction_9_KmS1, 2)*reaction_9_KmS2))/((reaction_9_KmI + species_3)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1)) + reaction_9_KmI*species_6*(1 + species_13/reaction_9_KmS2)*(reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/(pow(reaction_9_KmS1, 2)*(reaction_9_KmI + species_3)*pow((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1, 2));
            break;
        case 120:
            dwdp[20] = reaction_9_KmI*(species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/((reaction_9_KmI + species_3)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1));
            break;
        case 121:
            dwdp[20] = -reaction_9_KmI*(reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/(pow(reaction_9_KmI + species_3, 2)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1)) + (reaction_9_Vmax*species_13*species_6/(reaction_9_KmS1*reaction_9_KmS2) - reaction_9_Vmax*species_23*species_7/(reaction_9_Keq*reaction_9_KmS1*reaction_9_KmS2))/((reaction_9_KmI + species_3)*((1 + species_7/reaction_9_KmP1)*(1 + species_23/reaction_9_KmP2) + (1 + species_6/reaction_9_KmS1)*(1 + species_13/reaction_9_KmS2) - 1));
            break;
    }
}