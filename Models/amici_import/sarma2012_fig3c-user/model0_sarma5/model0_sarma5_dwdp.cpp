#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_sarma5(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*reaction_0_V1*species_0/(pow(reaction_0_K1, 2)*(1 + species_0/reaction_0_K1)*(1 + species_7/reaction_0_KI)) + 1.0*reaction_0_V1*pow(species_0, 2)/(pow(reaction_0_K1, 3)*pow(1 + species_0/reaction_0_K1, 2)*(1 + species_7/reaction_0_KI));
            break;
        case 1:
            dwdp[0] = 1.0*reaction_0_V1*species_0*species_7/(reaction_0_K1*pow(reaction_0_KI, 2)*(1 + species_0/reaction_0_K1)*pow(1 + species_7/reaction_0_KI, 2));
            break;
        case 2:
            dwdp[0] = 1.0*species_0/(reaction_0_K1*(1 + species_0/reaction_0_K1)*(1 + species_7/reaction_0_KI));
            break;
        case 3:
            dwdp[1] = 1.0*reaction_1_k3*species_1*species_2*species_7/(reaction_1_K3*reaction_1_Ka*(1 + species_7/reaction_1_Ka)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
            break;
        case 4:
            dwdp[1] = 1.0*reaction_1_k3*species_1*species_2*(species_2/pow(reaction_1_K3, 2) + species_3/pow(reaction_1_K3, 2))*(reaction_1_A*species_7/reaction_1_Ka + 1)/(reaction_1_K3*(1 + species_7/reaction_1_Ka)*pow(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3, 2)) - 1.0*reaction_1_k3*species_1*species_2*(reaction_1_A*species_7/reaction_1_Ka + 1)/(pow(reaction_1_K3, 2)*(1 + species_7/reaction_1_Ka)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
            break;
        case 5:
            dwdp[1] = -1.0*reaction_1_A*reaction_1_k3*species_1*species_2*species_7/(reaction_1_K3*pow(reaction_1_Ka, 2)*(1 + species_7/reaction_1_Ka)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3)) + 1.0*reaction_1_k3*species_1*species_2*species_7*(reaction_1_A*species_7/reaction_1_Ka + 1)/(reaction_1_K3*pow(reaction_1_Ka, 2)*pow(1 + species_7/reaction_1_Ka, 2)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
            break;
        case 6:
            dwdp[1] = 1.0*species_1*species_2*(reaction_1_A*species_7/reaction_1_Ka + 1)/(reaction_1_K3*(1 + species_7/reaction_1_Ka)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
            break;
        case 7:
            dwdp[2] = 1.0*reaction_2_k4*species_1*species_3*species_7/(reaction_2_K4*reaction_2_Ka*(1 + species_7/reaction_2_Ka)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
            break;
        case 8:
            dwdp[2] = 1.0*reaction_2_k4*species_1*species_3*(species_2/pow(reaction_2_K4, 2) + species_3/pow(reaction_2_K4, 2))*(reaction_2_A*species_7/reaction_2_Ka + 1)/(reaction_2_K4*(1 + species_7/reaction_2_Ka)*pow(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4, 2)) - 1.0*reaction_2_k4*species_1*species_3*(reaction_2_A*species_7/reaction_2_Ka + 1)/(pow(reaction_2_K4, 2)*(1 + species_7/reaction_2_Ka)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
            break;
        case 9:
            dwdp[2] = -1.0*reaction_2_A*reaction_2_k4*species_1*species_3*species_7/(reaction_2_K4*pow(reaction_2_Ka, 2)*(1 + species_7/reaction_2_Ka)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4)) + 1.0*reaction_2_k4*species_1*species_3*species_7*(reaction_2_A*species_7/reaction_2_Ka + 1)/(reaction_2_K4*pow(reaction_2_Ka, 2)*pow(1 + species_7/reaction_2_Ka, 2)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
            break;
        case 10:
            dwdp[2] = 1.0*species_1*species_3*(reaction_2_A*species_7/reaction_2_Ka + 1)/(reaction_2_K4*(1 + species_7/reaction_2_Ka)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
            break;
        case 11:
            dwdp[3] = 1.0*reaction_3_k7*species_4*species_5*(species_5/pow(reaction_3_K7, 2) + species_6/pow(reaction_3_K7, 2))/(reaction_3_K7*pow(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7, 2)) - 1.0*reaction_3_k7*species_4*species_5/(pow(reaction_3_K7, 2)*(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7));
            break;
        case 12:
            dwdp[3] = 1.0*species_4*species_5/(reaction_3_K7*(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7));
            break;
        case 13:
            dwdp[4] = 1.0*reaction_4_k8*species_4*species_6*(species_5/pow(reaction_4_K8, 2) + species_6/pow(reaction_4_K8, 2))/(reaction_4_K8*pow(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8, 2)) - 1.0*reaction_4_k8*species_4*species_6/(pow(reaction_4_K8, 2)*(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8));
            break;
        case 14:
            dwdp[4] = 1.0*species_4*species_6/(reaction_4_K8*(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8));
            break;
        case 15:
            dwdp[5] = -1.0*reaction_5_k2*species_1*species_8/(pow(reaction_5_K2, 2)*(1 + species_1/reaction_5_K2)) + 1.0*reaction_5_k2*pow(species_1, 2)*species_8/(pow(reaction_5_K2, 3)*pow(1 + species_1/reaction_5_K2, 2));
            break;
        case 16:
            dwdp[5] = 1.0*species_1*species_8/(reaction_5_K2*(1 + species_1/reaction_5_K2));
            break;
        case 17:
            dwdp[6] = 1.0*reaction_6_k5*species_4*species_9*(species_3/pow(reaction_6_K5, 2) + species_4/pow(reaction_6_K5, 2))/(reaction_6_K5*pow(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5, 2)) - 1.0*reaction_6_k5*species_4*species_9/(pow(reaction_6_K5, 2)*(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5));
            break;
        case 18:
            dwdp[6] = 1.0*species_4*species_9/(reaction_6_K5*(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5));
            break;
        case 19:
            dwdp[7] = 1.0*reaction_7_k6*species_3*species_9*(species_3/pow(reaction_7_K6, 2) + species_4/pow(reaction_7_K6, 2))/(reaction_7_K6*pow(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6, 2)) - 1.0*reaction_7_k6*species_3*species_9/(pow(reaction_7_K6, 2)*(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6));
            break;
        case 20:
            dwdp[7] = 1.0*species_3*species_9/(reaction_7_K6*(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6));
            break;
        case 21:
            dwdp[8] = 1.0*reaction_8_k9*species_10*species_7*(species_6/pow(reaction_8_K9, 2) + species_7/pow(reaction_8_K9, 2))/(reaction_8_K9*pow(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9, 2)) - 1.0*reaction_8_k9*species_10*species_7/(pow(reaction_8_K9, 2)*(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9));
            break;
        case 22:
            dwdp[8] = 1.0*species_10*species_7/(reaction_8_K9*(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9));
            break;
        case 23:
            dwdp[9] = 1.0*reaction_9_k10*species_10*species_6*(species_6/pow(reaction_9_K10, 2) + species_7/pow(reaction_9_K10, 2))/(reaction_9_K10*pow(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10, 2)) - 1.0*reaction_9_k10*species_10*species_6/(pow(reaction_9_K10, 2)*(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10));
            break;
        case 24:
            dwdp[9] = 1.0*species_10*species_6/(reaction_9_K10*(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10));
            break;
    }
}