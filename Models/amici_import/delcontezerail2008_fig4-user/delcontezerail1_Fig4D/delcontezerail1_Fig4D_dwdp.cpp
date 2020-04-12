#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_delcontezerail1_Fig4D(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0;
            break;
        case 1:
            dwdp[1] = -1.0*r5*reaction_1_ke*t*(-R5 + reaction_1_kg)*exp(reaction_1_kf*(-R5 + reaction_1_kg))/((t + 100)*pow(exp(reaction_1_kf*(-R5 + reaction_1_kg)) + 1, 2));
            break;
        case 2:
            dwdp[1] = -1.0*r5*reaction_1_ke*reaction_1_kf*t*exp(reaction_1_kf*(-R5 + reaction_1_kg))/((t + 100)*pow(exp(reaction_1_kf*(-R5 + reaction_1_kg)) + 1, 2));
            break;
        case 3:
            dwdp[1] = 1.0*r5*t/((t + 100)*(exp(reaction_1_kf*(-R5 + reaction_1_kg)) + 1));
            break;
        case 4:
            dwdp[2] = 1.0*r5;
            break;
        case 5:
            dwdp[3] = 1.0;
            break;
        case 6:
            dwdp[4] = -1.0*pow(R7, reaction_4_h)*r7*reaction_4_ke/pow(pow(R7, reaction_4_h) + reaction_4_kg, 2);
            break;
        case 7:
            dwdp[4] = -1.0*pow(R7, 2*reaction_4_h)*r7*reaction_4_ke*log(R7)/pow(pow(R7, reaction_4_h) + reaction_4_kg, 2) + 1.0*pow(R7, reaction_4_h)*r7*reaction_4_ke*log(R7)/(pow(R7, reaction_4_h) + reaction_4_kg);
            break;
        case 8:
            dwdp[4] = 1.0*pow(R7, reaction_4_h)*r7/(pow(R7, reaction_4_h) + reaction_4_kg);
            break;
        case 9:
            dwdp[5] = -1.0*r7*reaction_5_ke*(-R5 + reaction_5_kg)*exp(reaction_5_kf*(-R5 + reaction_5_kg))/pow(exp(reaction_5_kf*(-R5 + reaction_5_kg)) + 1, 2);
            break;
        case 10:
            dwdp[5] = -1.0*r7*reaction_5_ke*reaction_5_kf*exp(reaction_5_kf*(-R5 + reaction_5_kg))/pow(exp(reaction_5_kf*(-R5 + reaction_5_kg)) + 1, 2);
            break;
        case 11:
            dwdp[5] = 1.0*r7/(exp(reaction_5_kf*(-R5 + reaction_5_kg)) + 1);
            break;
        case 12:
            dwdp[6] = -1.0*R5*reaction_6_ke*(-R7 + reaction_6_kg)*exp(reaction_6_kf*(-R7 + reaction_6_kg))/pow(exp(reaction_6_kf*(-R7 + reaction_6_kg)) + 1, 2);
            break;
        case 13:
            dwdp[6] = -1.0*R5*reaction_6_ke*reaction_6_kf*exp(reaction_6_kf*(-R7 + reaction_6_kg))/pow(exp(reaction_6_kf*(-R7 + reaction_6_kg)) + 1, 2);
            break;
        case 14:
            dwdp[6] = 1.0*R5/(exp(reaction_6_kf*(-R7 + reaction_6_kg)) + 1);
            break;
        case 15:
            dwdp[7] = 1.0*r7;
            break;
        case 16:
            dwdp[8] = 1.0*R5;
            break;
        case 17:
            dwdp[9] = 1.0*R7;
            break;
    }
}