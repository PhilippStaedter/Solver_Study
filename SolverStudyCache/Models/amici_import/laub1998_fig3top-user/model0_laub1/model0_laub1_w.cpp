#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_laub1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*parameter_0*species_6;
    w[1] = 1.0*parameter_1*species_4;
    w[2] = 1.0*parameter_10*species_4;
    w[3] = 1.0*parameter_11*species_0;
    w[4] = 1.0*parameter_12*species_0;
    w[5] = 1.0*parameter_13*species_2*species_5;
    w[6] = 1.0*parameter_2*species_1;
    w[7] = 1.0*parameter_3*species_2;
    w[8] = 1.0*parameter_4*species_5;
    w[9] = 1.0*parameter_5*species_2*species_6;
    w[10] = 1.0*parameter_6;
    w[11] = 1.0*parameter_7*species_3*species_6;
    w[12] = 1.0*parameter_8*species_4;
    w[13] = 1.0*parameter_9*species_1*species_3;
}