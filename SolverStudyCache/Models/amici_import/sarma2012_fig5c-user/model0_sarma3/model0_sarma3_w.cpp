#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_sarma3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*reaction_1_k1*species_1*species_11/(parameter_1 + species_1);
    w[1] = 1.0*reaction_10_k10b*species_10*species_7/(parameter_10*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    w[2] = 1.0*reaction_2_k2a*species_2*species_9/(parameter_2*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11));
    w[3] = 1.0*reaction_3_k3*species_2*species_3/(parameter_3*(1 + species_4/parameter_4 + species_3/parameter_3));
    w[4] = 1.0*reaction_4_k4*species_2*species_4/(parameter_4*(1 + species_4/parameter_4 + species_3/parameter_3));
    w[5] = 1.0*reaction_5_k5a*species_5*species_9/(parameter_5*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) + 1.0*reaction_5_k5b*species_10*species_5/(parameter_13*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    w[6] = 1.0*reaction_6_k6a*species_4*species_9/(parameter_6*(1 + species_4/parameter_6 + species_5/parameter_5 + species_2/parameter_2 + species_1/parameter_11 + species_3/parameter_11)) + 1.0*reaction_6_k6b*species_10*species_4/(parameter_14*(1 + species_8/parameter_9 + species_4/parameter_14 + species_5/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
    w[7] = 1.0*reaction_7_k7*species_5*species_6/(parameter_7*(1 + species_7/parameter_8 + species_6/parameter_7));
    w[8] = 1.0*reaction_8_k7*species_5*species_7/(parameter_8*(1 + species_7/parameter_8 + species_6/parameter_7));
    w[9] = 1.0*reaction_9_k9b*species_10*species_8/(parameter_9*(1 + species_8/parameter_9 + species_4/parameter_14 + species_9/parameter_13 + species_3/parameter_12 + species_6/parameter_12 + species_7/parameter_10));
}