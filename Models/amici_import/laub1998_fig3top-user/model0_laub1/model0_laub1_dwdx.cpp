#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_laub1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*parameter_11;
    dwdx[1] = 1.0*parameter_12;
    dwdx[2] = 1.0*parameter_2;
    dwdx[3] = 1.0*parameter_9*species_3;
    dwdx[4] = 1.0*parameter_13*species_5;
    dwdx[5] = 1.0*parameter_3;
    dwdx[6] = 1.0*parameter_5*species_6;
    dwdx[7] = 1.0*parameter_7*species_6;
    dwdx[8] = 1.0*parameter_9*species_1;
    dwdx[9] = 1.0*parameter_1;
    dwdx[10] = 1.0*parameter_10;
    dwdx[11] = 1.0*parameter_8;
    dwdx[12] = 1.0*parameter_13*species_2;
    dwdx[13] = 1.0*parameter_4;
    dwdx[14] = 1.0*parameter_0;
    dwdx[15] = 1.0*parameter_5*species_2;
    dwdx[16] = 1.0*parameter_7*species_3;
}