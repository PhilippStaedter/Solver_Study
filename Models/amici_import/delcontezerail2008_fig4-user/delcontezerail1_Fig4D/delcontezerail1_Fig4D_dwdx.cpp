#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_delcontezerail1_Fig4D(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*r5*reaction_1_ke*reaction_1_kf*t*exp(reaction_1_kf*(-R5 + reaction_1_kg))/((t + 100)*pow(exp(reaction_1_kf*(-R5 + reaction_1_kg)) + 1, 2));
    dwdx[1] = 1.0*r7*reaction_5_ke*reaction_5_kf*exp(reaction_5_kf*(-R5 + reaction_5_kg))/pow(exp(reaction_5_kf*(-R5 + reaction_5_kg)) + 1, 2);
    dwdx[2] = 1.0*reaction_6_ke/(exp(reaction_6_kf*(-R7 + reaction_6_kg)) + 1);
    dwdx[3] = 1.0*reaction_8_kh;
    dwdx[4] = -1.0*pow(R7, 2*reaction_4_h)*r7*reaction_4_h*reaction_4_ke/(R7*pow(pow(R7, reaction_4_h) + reaction_4_kg, 2)) + 1.0*pow(R7, reaction_4_h)*r7*reaction_4_h*reaction_4_ke/(R7*(pow(R7, reaction_4_h) + reaction_4_kg));
    dwdx[5] = 1.0*R5*reaction_6_ke*reaction_6_kf*exp(reaction_6_kf*(-R7 + reaction_6_kg))/pow(exp(reaction_6_kf*(-R7 + reaction_6_kg)) + 1, 2);
    dwdx[6] = 1.0*reaction_9_kh;
    dwdx[7] = 1.0*reaction_1_ke*t/((t + 100)*(exp(reaction_1_kf*(-R5 + reaction_1_kg)) + 1));
    dwdx[8] = 1.0*reaction_2_kminus1;
    dwdx[9] = 1.0*pow(R7, reaction_4_h)*reaction_4_ke/(pow(R7, reaction_4_h) + reaction_4_kg);
    dwdx[10] = 1.0*reaction_5_ke/(exp(reaction_5_kf*(-R5 + reaction_5_kg)) + 1);
    dwdx[11] = 1.0*reaction_7_kminus1;
}