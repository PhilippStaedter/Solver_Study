#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kim5(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*reaction_1_k1;
    dwdx[1] = -1.0*parameter_5*reaction_8_V*pow(species_1, 2*parameter_5)/(species_1*pow(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5), 2)) + 1.0*parameter_5*reaction_8_V*pow(species_1, parameter_5)/(species_1*(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5)));
    dwdx[2] = 1.0*reaction_2_k1;
    dwdx[3] = -1.0*parameter_4*reaction_7_V*pow(species_2, parameter_4)/(species_2*pow(pow(reaction_7_Shalve, parameter_4) + pow(species_2, parameter_4), 2));
    dwdx[4] = 1.0*reaction_3_k1;
    dwdx[5] = 1.0*parameter_1;
    dwdx[6] = 1.0*reaction_4_k1;
    dwdx[7] = 1.0*parameter_2;
    dwdx[8] = 1.0*parameter_6;
}