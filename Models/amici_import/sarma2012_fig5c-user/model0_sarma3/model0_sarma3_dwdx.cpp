#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_sarma3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*reaction_1_k1*species_1*species_11/pow(parameter_1 + species_1, 2) + 1.0*reaction_1_k1*species_11/(parameter_1 + species_1);
    dwdx[1] = -1.0*reaction_2_k2a*species_2*species_9/(parameter_11*parameter_2*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[2] = -1.0*reaction_5_k5a*species_5*species_9/(parameter_11*parameter_5*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[3] = -1.0*reaction_6_k6a*species_4*species_9/(parameter_11*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[4] = 1.0*reaction_10_k10b*species_7/(parameter_10*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    dwdx[5] = 1.0*reaction_5_k5b*species_5/(parameter_13*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    dwdx[6] = 1.0*reaction_6_k6b*species_4/(parameter_14*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    dwdx[7] = 1.0*reaction_9_k9b*species_8/(parameter_9*(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    dwdx[8] = 1.0*reaction_1_k1*species_1/(parameter_1 + species_1);
    dwdx[9] = 1.0*reaction_2_k2a*species_9/(parameter_2*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) - 1.0*reaction_2_k2a*species_2*species_9/(pow(parameter_2, 2)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[10] = 1.0*reaction_3_k3*species_3/(parameter_3*(1 + species_4/parameter_4 + species_3/parameter_3));
    dwdx[11] = 1.0*reaction_4_k4*species_4/(parameter_4*(1 + species_4/parameter_4 + species_3/parameter_3));
    dwdx[12] = -1.0*reaction_5_k5a*species_5*species_9/(parameter_2*parameter_5*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[13] = -1.0*reaction_6_k6a*species_4*species_9/(parameter_2*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[14] = -1.0*reaction_10_k10b*species_10*species_7/(parameter_10*parameter_12*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[15] = -1.0*reaction_2_k2a*species_2*species_9/(parameter_11*parameter_2*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[16] = 1.0*reaction_3_k3*species_2/(parameter_3*(1 + species_4/parameter_4 + species_3/parameter_3)) - 1.0*reaction_3_k3*species_2*species_3/(pow(parameter_3, 2)*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
    dwdx[17] = -1.0*reaction_4_k4*species_2*species_4/(parameter_3*parameter_4*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
    dwdx[18] = -1.0*reaction_5_k5b*species_10*species_5/(parameter_12*parameter_13*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2)) - 1.0*reaction_5_k5a*species_5*species_9/(parameter_11*parameter_5*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[19] = -1.0*reaction_6_k6b*species_10*species_4/(parameter_12*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2)) - 1.0*reaction_6_k6a*species_4*species_9/(parameter_11*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[20] = -1.0*reaction_9_k9b*species_10*species_8/(parameter_12*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[21] = -1.0*reaction_10_k10b*species_10*species_7/(parameter_10*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[22] = -1.0*reaction_2_k2a*species_2*species_9/(parameter_2*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[23] = -1.0*reaction_3_k3*species_2*species_3/(parameter_3*parameter_4*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
    dwdx[24] = 1.0*reaction_4_k4*species_2/(parameter_4*(1 + species_4/parameter_4 + species_3/parameter_3)) - 1.0*reaction_4_k4*species_2*species_4/(pow(parameter_4, 2)*pow(1 + species_4/parameter_4 + species_3/parameter_3, 2));
    dwdx[25] = -1.0*reaction_5_k5a*species_5*species_9/(parameter_5*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2)) - 1.0*reaction_5_k5b*species_10*species_5/(parameter_13*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[26] = 1.0*reaction_6_k6a*species_9/(parameter_6*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) - 1.0*reaction_6_k6a*species_4*species_9/(pow(parameter_6, 2)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2)) + 1.0*reaction_6_k6b*species_10/(parameter_14*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) - 1.0*reaction_6_k6b*species_10*species_4/(pow(parameter_14, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[27] = -1.0*reaction_9_k9b*species_10*species_8/(parameter_14*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[28] = -1.0*reaction_10_k10b*species_10*species_7/(parameter_10*parameter_13*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[29] = -1.0*reaction_2_k2a*species_2*species_9/(parameter_2*parameter_5*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2));
    dwdx[30] = 1.0*reaction_5_k5a*species_9/(parameter_5*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) - 1.0*reaction_5_k5a*species_5*species_9/(pow(parameter_5, 2)*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2)) + 1.0*reaction_5_k5b*species_10/(parameter_13*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) - 1.0*reaction_5_k5b*species_10*species_5/(pow(parameter_13, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[31] = -1.0*reaction_6_k6a*species_4*species_9/(parameter_5*parameter_6*pow(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11, 2)) - 1.0*reaction_6_k6b*species_10*species_4/(parameter_13*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[32] = 1.0*reaction_7_k7*species_6/(parameter_7*(1 + species_7/parameter_8 + species_6/parameter_7));
    dwdx[33] = 1.0*reaction_8_k7*species_7/(parameter_8*(1 + species_7/parameter_8 + species_6/parameter_7));
    dwdx[34] = -1.0*reaction_10_k10b*species_10*species_7/(parameter_10*parameter_12*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[35] = -1.0*reaction_5_k5b*species_10*species_5/(parameter_12*parameter_13*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[36] = -1.0*reaction_6_k6b*species_10*species_4/(parameter_12*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[37] = 1.0*reaction_7_k7*species_5/(parameter_7*(1 + species_7/parameter_8 + species_6/parameter_7)) - 1.0*reaction_7_k7*species_5*species_6/(pow(parameter_7, 2)*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
    dwdx[38] = -1.0*reaction_8_k7*species_5*species_7/(parameter_7*parameter_8*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
    dwdx[39] = -1.0*reaction_9_k9b*species_10*species_8/(parameter_12*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[40] = 1.0*reaction_10_k10b*species_10/(parameter_10*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) - 1.0*reaction_10_k10b*species_10*species_7/(pow(parameter_10, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[41] = -1.0*reaction_5_k5b*species_10*species_5/(parameter_10*parameter_13*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[42] = -1.0*reaction_6_k6b*species_10*species_4/(parameter_10*parameter_14*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[43] = -1.0*reaction_7_k7*species_5*species_6/(parameter_7*parameter_8*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
    dwdx[44] = 1.0*reaction_8_k7*species_5/(parameter_8*(1 + species_7/parameter_8 + species_6/parameter_7)) - 1.0*reaction_8_k7*species_5*species_7/(pow(parameter_8, 2)*pow(1 + species_7/parameter_8 + species_6/parameter_7, 2));
    dwdx[45] = -1.0*reaction_9_k9b*species_10*species_8/(parameter_10*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[46] = -1.0*reaction_10_k10b*species_10*species_7/(parameter_10*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[47] = -1.0*reaction_5_k5b*species_10*species_5/(parameter_13*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[48] = -1.0*reaction_6_k6b*species_10*species_4/(parameter_14*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[49] = 1.0*reaction_9_k9b*species_10/(parameter_9*(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10)) - 1.0*reaction_9_k9b*species_10*species_8/(pow(parameter_9, 2)*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
    dwdx[50] = 1.0*reaction_2_k2a*species_2/(parameter_2*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
    dwdx[51] = 1.0*reaction_5_k5a*species_5/(parameter_5*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
    dwdx[52] = 1.0*reaction_6_k6a*species_4/(parameter_6*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
    dwdx[53] = -1.0*reaction_9_k9b*species_10*species_8/(parameter_13*parameter_9*pow(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10, 2));
}