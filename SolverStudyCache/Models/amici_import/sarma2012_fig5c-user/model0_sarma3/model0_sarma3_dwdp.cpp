#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_sarma3(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*reaction_1_k1*species_1*species_11/pow(parameter_1 + species_1, 2);
            break;
        case 1:
            dwdp[1] = -1.0*reaction_10_k10b*species_10*species_7/(pow(parameter_10, 2)*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) + 1.0*reaction_10_k10b*species_10*pow(species_7, 2)/(pow(parameter_10, 3)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[5] = 1.0*reaction_5_k5b*species_10*species_5*species_7/(pow(parameter_10, 2)*parameter_13*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[6] = 1.0*reaction_6_k6b*species_10*species_4*species_7/(pow(parameter_10, 2)*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[9] = 1.0*reaction_9_k9b*species_10*species_7*species_8/(pow(parameter_10, 2)*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            break;
        case 2:
            dwdp[2] = 1.0*reaction_2_k2a*species_2*species_9*(species_1/pow(parameter_11, 2) + species_3/pow(parameter_11, 2))/(parameter_2*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[5] = 1.0*reaction_5_k5a*species_5*species_9*(species_1/pow(parameter_11, 2) + species_3/pow(parameter_11, 2))/(parameter_5*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[6] = 1.0*reaction_6_k6a*species_4*species_9*(species_1/pow(parameter_11, 2) + species_3/pow(parameter_11, 2))/(parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            break;
        case 3:
            dwdp[1] = 1.0*reaction_10_k10b*species_10*species_7*(species_3/pow(parameter_12, 2) + species_6/pow(parameter_12, 2))/(parameter_10*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[5] = 1.0*reaction_5_k5b*species_10*species_5*(species_3/pow(parameter_12, 2) + species_6/pow(parameter_12, 2))/(parameter_13*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[6] = 1.0*reaction_6_k6b*species_10*species_4*(species_3/pow(parameter_12, 2) + species_6/pow(parameter_12, 2))/(parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[9] = 1.0*reaction_9_k9b*species_10*species_8*(species_3/pow(parameter_12, 2) + species_6/pow(parameter_12, 2))/(parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            break;
        case 4:
            dwdp[1] = 1.0*reaction_10_k10b*species_10*species_5*species_7/(parameter_10*pow(parameter_13, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[5] = -1.0*reaction_5_k5b*species_10*species_5/(pow(parameter_13, 2)*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) + 1.0*reaction_5_k5b*species_10*pow(species_5, 2)/(pow(parameter_13, 3)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[6] = 1.0*reaction_6_k6b*species_10*species_4*species_5/(pow(parameter_13, 2)*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[9] = 1.0*reaction_9_k9b*species_10*species_8*species_9/(pow(parameter_13, 2)*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            break;
        case 5:
            dwdp[1] = 1.0*reaction_10_k10b*species_10*species_4*species_7/(parameter_10*pow(parameter_14, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[5] = 1.0*reaction_5_k5b*species_10*species_4*species_5/(parameter_13*pow(parameter_14, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[6] = -1.0*reaction_6_k6b*species_10*species_4/(pow(parameter_14, 2)*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) + 1.0*reaction_6_k6b*species_10*pow(species_4, 2)/(pow(parameter_14, 3)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[9] = 1.0*reaction_9_k9b*species_10*species_4*species_8/(pow(parameter_14, 2)*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            break;
        case 6:
            dwdp[2] = -1.0*reaction_2_k2a*species_2*species_9/(pow(parameter_2, 2)*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) + 1.0*reaction_2_k2a*pow(species_2, 2)*species_9/(pow(parameter_2, 3)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[5] = 1.0*reaction_5_k5a*species_2*species_5*species_9/(pow(parameter_2, 2)*parameter_5*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[6] = 1.0*reaction_6_k6a*species_2*species_4*species_9/(pow(parameter_2, 2)*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            break;
        case 7:
            dwdp[3] = -1.0*reaction_3_k3*species_2*species_3/(pow(parameter_3, 2)*(1 + species_4/parameter_4 + species_3/parameter_3)) + 1.0*reaction_3_k3*species_2*pow(species_3, 2)/(pow(parameter_3, 3)*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
            dwdp[4] = 1.0*reaction_4_k4*species_2*species_3*species_4/(pow(parameter_3, 2)*parameter_4*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
            break;
        case 8:
            dwdp[3] = 1.0*reaction_3_k3*species_2*species_3*species_4/(parameter_3*pow(parameter_4, 2)*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
            dwdp[4] = -1.0*reaction_4_k4*species_2*species_4/(pow(parameter_4, 2)*(1 + species_4/parameter_4 + species_3/parameter_3)) + 1.0*reaction_4_k4*species_2*pow(species_4, 2)/(pow(parameter_4, 3)*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
            break;
        case 9:
            dwdp[2] = 1.0*reaction_2_k2a*species_2*species_5*species_9/(parameter_2*pow(parameter_5, 2)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[5] = -1.0*reaction_5_k5a*species_5*species_9/(pow(parameter_5, 2)*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) + 1.0*reaction_5_k5a*pow(species_5, 2)*species_9/(pow(parameter_5, 3)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[6] = 1.0*reaction_6_k6a*species_4*species_5*species_9/(pow(parameter_5, 2)*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            break;
        case 10:
            dwdp[2] = 1.0*reaction_2_k2a*species_2*species_4*species_9/(parameter_2*pow(parameter_6, 2)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[5] = 1.0*reaction_5_k5a*species_4*species_5*species_9/(parameter_5*pow(parameter_6, 2)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            dwdp[6] = -1.0*reaction_6_k6a*species_4*species_9/(pow(parameter_6, 2)*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) + 1.0*reaction_6_k6a*pow(species_4, 2)*species_9/(pow(parameter_6, 3)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
            break;
        case 11:
            dwdp[7] = -1.0*reaction_7_k7*species_5*species_6/(pow(parameter_7, 2)*(1 + species_7/parameter_8 + species_6/parameter_7)) + 1.0*reaction_7_k7*species_5*pow(species_6, 2)/(pow(parameter_7, 3)*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
            dwdp[8] = 1.0*reaction_8_k7*species_5*species_6*species_7/(pow(parameter_7, 2)*parameter_8*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
            break;
        case 12:
            dwdp[7] = 1.0*reaction_7_k7*species_5*species_6*species_7/(parameter_7*pow(parameter_8, 2)*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
            dwdp[8] = -1.0*reaction_8_k7*species_5*species_7/(pow(parameter_8, 2)*(1 + species_7/parameter_8 + species_6/parameter_7)) + 1.0*reaction_8_k7*species_5*pow(species_7, 2)/(pow(parameter_8, 3)*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
            break;
        case 13:
            dwdp[1] = 1.0*reaction_10_k10b*species_10*species_7*species_8/(parameter_10*pow(parameter_9, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[5] = 1.0*reaction_5_k5b*species_10*species_5*species_8/(parameter_13*pow(parameter_9, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[6] = 1.0*reaction_6_k6b*species_10*species_4*species_8/(parameter_14*pow(parameter_9, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            dwdp[9] = -1.0*reaction_9_k9b*species_10*species_8/(pow(parameter_9, 2)*(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) + 1.0*reaction_9_k9b*species_10*pow(species_8, 2)/(pow(parameter_9, 3)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
            break;
        case 14:
            dwdp[0] = 1.0*species_1*species_11/(parameter_1 + species_1);
            break;
        case 15:
            dwdp[1] = 1.0*species_10*species_7/(parameter_10*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
            break;
        case 16:
            dwdp[2] = 1.0*species_2*species_9/(parameter_2*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
            break;
        case 17:
            dwdp[3] = 1.0*species_2*species_3/(parameter_3*(1 + species_4/parameter_4 + species_3/parameter_3));
            break;
        case 18:
            dwdp[4] = 1.0*species_2*species_4/(parameter_4*(1 + species_4/parameter_4 + species_3/parameter_3));
            break;
        case 19:
            dwdp[5] = 1.0*species_5*species_9/(parameter_5*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
            break;
        case 20:
            dwdp[5] = 1.0*species_10*species_5/(parameter_13*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
            break;
        case 21:
            dwdp[6] = 1.0*species_4*species_9/(parameter_6*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
            break;
        case 22:
            dwdp[6] = 1.0*species_10*species_4/(parameter_14*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
            break;
        case 23:
            dwdp[7] = 1.0*species_5*species_6/(parameter_7*(1 + species_7/parameter_8 + species_6/parameter_7));
            break;
        case 24:
            dwdp[8] = 1.0*species_5*species_7/(parameter_8*(1 + species_7/parameter_8 + species_6/parameter_7));
            break;
        case 25:
            dwdp[9] = 1.0*species_10*species_8/(parameter_9*(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
            break;
    }
}