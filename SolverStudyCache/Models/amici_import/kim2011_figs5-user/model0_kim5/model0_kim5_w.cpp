#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kim5(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*reaction_1_k1*species_1;
    w[1] = 1.0*reaction_2_k1*species_2;
    w[2] = 1.0*reaction_3_k1*species_3;
    w[3] = 1.0*reaction_4_k1*species_4;
    w[4] = 1.0*parameter_1*species_3;
    w[5] = 1.0*parameter_2*species_4;
    w[6] = 1.0*reaction_7_V/(pow(reaction_7_Shalve, parameter_4) + pow(species_2, parameter_4));
    w[7] = 1.0*reaction_8_V*pow(species_1, parameter_5)/(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5));
    w[8] = 1.0*parameter_6*species_4;
}