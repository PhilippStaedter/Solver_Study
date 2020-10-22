#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_alexander1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*A*b1;
    w[1] = 1.0*A*R*sigma1;
    w[2] = 1.0*G*v;
    w[3] = 1.0*G*f*v;
    w[4] = 1.0*E*gamma;
    w[5] = 1.0*A*beta;
    w[6] = 1.0*A*E*pi1;
    w[7] = 1.0*A*lambdaE;
    w[8] = 1.0*A*muA;
    w[9] = 1.0*R*muR;
    w[10] = 1.0*E*muE;
    w[11] = 1.0*G*muG;
}