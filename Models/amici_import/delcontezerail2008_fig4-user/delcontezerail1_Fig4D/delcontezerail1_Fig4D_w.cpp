#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_delcontezerail1_Fig4D(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*reaction_0_K1;
    w[1] = 1.0*r5*reaction_1_ke*t/((t + 100)*(exp(reaction_1_kf*(-R5 + reaction_1_kg)) + 1));
    w[2] = 1.0*r5*reaction_2_kminus1;
    w[3] = 1.0*reaction_3_K1;
    w[4] = 1.0*pow(R7, reaction_4_h)*r7*reaction_4_ke/(pow(R7, reaction_4_h) + reaction_4_kg);
    w[5] = 1.0*r7*reaction_5_ke/(exp(reaction_5_kf*(-R5 + reaction_5_kg)) + 1);
    w[6] = 1.0*R5*reaction_6_ke/(exp(reaction_6_kf*(-R7 + reaction_6_kg)) + 1);
    w[7] = 1.0*r7*reaction_7_kminus1;
    w[8] = 1.0*R5*reaction_8_kh;
    w[9] = 1.0*R7*reaction_9_kh;
}