#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model1_alexander2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*b1;
    dwdx[1] = 1.0*R*sigma1;
    dwdx[2] = 1.0*beta;
    dwdx[3] = 1.0*E*pi1;
    dwdx[4] = 1.0*lambdaE;
    dwdx[5] = 1.0*muA;
    dwdx[6] = 1.0*gamma;
    dwdx[7] = 1.0*A*pi1;
    dwdx[8] = 1.0*muE;
    dwdx[9] = -1.0*G*v_max/pow(1.0*G + amici_k, 2) + 1.0*v_max/(1.0*G + amici_k);
    dwdx[10] = -1.0*G*f*v_max/pow(1.0*G + amici_k, 2) + 1.0*f*v_max/(1.0*G + amici_k);
    dwdx[11] = 1.0*muG;
    dwdx[12] = 1.0*A*sigma1;
    dwdx[13] = 1.0*muR;
}