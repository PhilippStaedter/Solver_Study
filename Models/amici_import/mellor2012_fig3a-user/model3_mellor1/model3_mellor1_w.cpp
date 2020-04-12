#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model3_mellor1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1000.0*parameter_2*species_1/(parameter_1 + species_1);
    w[1] = 1000.0*parameter_4*species_1/(parameter_3 + species_1);
    w[2] = 1000.0*parameter_6*species_1/(parameter_5 + species_1);
    w[3] = 1000.0*parameter_8*species_7/(parameter_7 + species_7);
    w[4] = 135.0*parameter_8*species_8/(parameter_7 + species_8);
}