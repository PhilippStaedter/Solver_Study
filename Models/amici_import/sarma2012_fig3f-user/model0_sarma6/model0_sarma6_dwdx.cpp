#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_sarma6(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*reaction_0_V1*(reaction_0_A*species_7/reaction_0_Ka + 1)/(reaction_0_K1*(1 + species_0/reaction_0_K1)*(1 + species_7/reaction_0_Ka)) - 1.0*reaction_0_V1*species_0*(reaction_0_A*species_7/reaction_0_Ka + 1)/(pow(reaction_0_K1, 2)*pow(1 + species_0/reaction_0_K1, 2)*(1 + species_7/reaction_0_Ka));
    dwdx[1] = 1.0*reaction_1_k3*species_2/(reaction_1_K3*(1 + species_7/reaction_1_KI)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
    dwdx[2] = 1.0*reaction_2_k4*species_3/(reaction_2_K4*(1 + species_7/reaction_2_KI)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
    dwdx[3] = 1.0*reaction_5_k2*species_8/(reaction_5_K2*(1 + species_1/reaction_5_K2)) - 1.0*reaction_5_k2*species_1*species_8/(pow(reaction_5_K2, 2)*pow(1 + species_1/reaction_5_K2, 2));
    dwdx[4] = 1.0*reaction_8_k9*species_7/(reaction_8_K9*(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9));
    dwdx[5] = 1.0*reaction_9_k10*species_6/(reaction_9_K10*(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10));
    dwdx[6] = 1.0*reaction_1_k3*species_1/(reaction_1_K3*(1 + species_7/reaction_1_KI)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3)) - 1.0*reaction_1_k3*species_1*species_2/(pow(reaction_1_K3, 2)*(1 + species_7/reaction_1_KI)*pow(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3, 2));
    dwdx[7] = -1.0*reaction_2_k4*species_1*species_3/(pow(reaction_2_K4, 2)*(1 + species_7/reaction_2_KI)*pow(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4, 2));
    dwdx[8] = -1.0*reaction_1_k3*species_1*species_2/(pow(reaction_1_K3, 2)*(1 + species_7/reaction_1_KI)*pow(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3, 2));
    dwdx[9] = 1.0*reaction_2_k4*species_1/(reaction_2_K4*(1 + species_7/reaction_2_KI)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4)) - 1.0*reaction_2_k4*species_1*species_3/(pow(reaction_2_K4, 2)*(1 + species_7/reaction_2_KI)*pow(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4, 2));
    dwdx[10] = -1.0*reaction_6_k5*species_4*species_9/(pow(reaction_6_K5, 2)*pow(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5, 2));
    dwdx[11] = 1.0*reaction_7_k6*species_9/(reaction_7_K6*(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6)) - 1.0*reaction_7_k6*species_3*species_9/(pow(reaction_7_K6, 2)*pow(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6, 2));
    dwdx[12] = 1.0*reaction_3_k7*species_5/(reaction_3_K7*(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7));
    dwdx[13] = 1.0*reaction_4_k8*species_6/(reaction_4_K8*(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8));
    dwdx[14] = 1.0*reaction_6_k5*species_9/(reaction_6_K5*(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5)) - 1.0*reaction_6_k5*species_4*species_9/(pow(reaction_6_K5, 2)*pow(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5, 2));
    dwdx[15] = -1.0*reaction_7_k6*species_3*species_9/(pow(reaction_7_K6, 2)*pow(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6, 2));
    dwdx[16] = 1.0*reaction_3_k7*species_4/(reaction_3_K7*(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7)) - 1.0*reaction_3_k7*species_4*species_5/(pow(reaction_3_K7, 2)*pow(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7, 2));
    dwdx[17] = -1.0*reaction_4_k8*species_4*species_6/(pow(reaction_4_K8, 2)*pow(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8, 2));
    dwdx[18] = -1.0*reaction_3_k7*species_4*species_5/(pow(reaction_3_K7, 2)*pow(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7, 2));
    dwdx[19] = 1.0*reaction_4_k8*species_4/(reaction_4_K8*(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8)) - 1.0*reaction_4_k8*species_4*species_6/(pow(reaction_4_K8, 2)*pow(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8, 2));
    dwdx[20] = -1.0*reaction_8_k9*species_10*species_7/(pow(reaction_8_K9, 2)*pow(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9, 2));
    dwdx[21] = 1.0*reaction_9_k10*species_10/(reaction_9_K10*(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10)) - 1.0*reaction_9_k10*species_10*species_6/(pow(reaction_9_K10, 2)*pow(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10, 2));
    dwdx[22] = 1.0*reaction_0_A*reaction_0_V1*species_0/(reaction_0_K1*reaction_0_Ka*(1 + species_0/reaction_0_K1)*(1 + species_7/reaction_0_Ka)) - 1.0*reaction_0_V1*species_0*(reaction_0_A*species_7/reaction_0_Ka + 1)/(reaction_0_K1*reaction_0_Ka*(1 + species_0/reaction_0_K1)*pow(1 + species_7/reaction_0_Ka, 2));
    dwdx[23] = -1.0*reaction_1_k3*species_1*species_2/(reaction_1_K3*reaction_1_KI*pow(1 + species_7/reaction_1_KI, 2)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
    dwdx[24] = -1.0*reaction_2_k4*species_1*species_3/(reaction_2_K4*reaction_2_KI*pow(1 + species_7/reaction_2_KI, 2)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
    dwdx[25] = 1.0*reaction_8_k9*species_10/(reaction_8_K9*(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9)) - 1.0*reaction_8_k9*species_10*species_7/(pow(reaction_8_K9, 2)*pow(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9, 2));
    dwdx[26] = -1.0*reaction_9_k10*species_10*species_6/(pow(reaction_9_K10, 2)*pow(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10, 2));
    dwdx[27] = 1.0*reaction_5_k2*species_1/(reaction_5_K2*(1 + species_1/reaction_5_K2));
    dwdx[28] = 1.0*reaction_6_k5*species_4/(reaction_6_K5*(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5));
    dwdx[29] = 1.0*reaction_7_k6*species_3/(reaction_7_K6*(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6));
}