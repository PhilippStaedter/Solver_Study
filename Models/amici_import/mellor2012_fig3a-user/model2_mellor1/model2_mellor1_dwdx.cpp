#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_mellor1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1000.0*parameter_2*species_1/pow(parameter_1 + species_1, 2) + 1000.0*parameter_2/(parameter_1 + species_1);
    dwdx[1] = -1000.0*parameter_4*species_1/pow(parameter_3 + species_1, 2) + 1000.0*parameter_4/(parameter_3 + species_1);
    dwdx[2] = -1000.0*parameter_6*species_1/pow(parameter_5 + species_1, 2) + 1000.0*parameter_6/(parameter_5 + species_1);
    dwdx[3] = -1000.0*parameter_8*species_7/pow(parameter_7 + species_7, 2) + 1000.0*parameter_8/(parameter_7 + species_7);
    dwdx[4] = -135.0*parameter_8*species_8/pow(parameter_7 + species_8, 2) + 135.0*parameter_8/(parameter_7 + species_8);
}