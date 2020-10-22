#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_piedrafita1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*S*STU*k1 - 1.0*STUS*k1r;
    w[1] = -1.0*STU*SU*k10r + 1.0*STUSU*k10;
    w[2] = 1.0*ST*k4;
    w[3] = 1.0*STUS*T*k2 - 1.0*STUST*k2r;
    w[4] = -1.0*ST*STU*k3r + 1.0*STUST*k3;
    w[5] = 1.0*STU*k4;
    w[6] = 1.0*ST*SU*k5 - 1.0*SUST*k5r;
    w[7] = 1.0*SUST*U*k6 - 1.0*SUSTU*k6r;
    w[8] = -1.0*STU*SU*k7r + 1.0*SUSTU*k7;
    w[9] = 1.0*SU*k4;
    w[10] = 1.0*STUS*U*k9 - 1.0*STUSU*k9r;
}