#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_sarma7(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*reaction_0_V1*species_0*(reaction_0_A*species_7/reaction_0_Ka + 1)/(reaction_0_K1*(1 + species_0/reaction_0_K1)*(1 + species_7/reaction_0_Ka));
    w[1] = 1.0*reaction_1_k3*species_1*species_2/(reaction_1_K3*(1 + species_7/reaction_1_KI)*(1 + species_2/reaction_1_K3 + species_3/reaction_1_K3));
    w[2] = 1.0*reaction_10_V11*species_11*(reaction_10_A*species_7/reaction_10_Ka + 1)/(reaction_10_K11*(1 + species_11/reaction_10_K11)*(1 + species_7/reaction_10_Ka));
    w[3] = 1.0*reaction_11_V12*species_12/(reaction_11_K12 + species_12);
    w[4] = 1.0*reaction_2_k4*species_1*species_3/(reaction_2_K4*(1 + species_7/reaction_2_KI)*(1 + species_2/reaction_2_K4 + species_3/reaction_2_K4));
    w[5] = 1.0*reaction_3_k7*species_4*species_5/(reaction_3_K7*(1 + species_5/reaction_3_K7 + species_6/reaction_3_K7));
    w[6] = 1.0*reaction_4_k8*species_4*species_6/(reaction_4_K8*(1 + species_5/reaction_4_K8 + species_6/reaction_4_K8));
    w[7] = 1.0*reaction_5_k2*species_1*species_8/(reaction_5_K2*(1 + species_1/reaction_5_K2));
    w[8] = 1.0*reaction_6_k5*species_4*species_9/(reaction_6_K5*(1 + species_3/reaction_6_K5 + species_4/reaction_6_K5));
    w[9] = 1.0*reaction_7_k6*species_3*species_9/(reaction_7_K6*(1 + species_3/reaction_7_K6 + species_4/reaction_7_K6));
    w[10] = 1.0*reaction_8_k9*species_10*species_7/(reaction_8_K9*(1 + species_6/reaction_8_K9 + species_7/reaction_8_K9));
    w[11] = 1.0*reaction_9_k10*species_10*species_6/(reaction_9_K10*(1 + species_6/reaction_9_K10 + species_7/reaction_9_K10));
}