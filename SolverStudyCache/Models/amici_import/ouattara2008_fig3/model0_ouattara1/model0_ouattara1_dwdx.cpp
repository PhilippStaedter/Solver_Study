#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_ouattara1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = delta;
    dwdx[1] = V*beta;
    dwdx[2] = mu;
    dwdx[3] = amici_k;
    dwdx[4] = T*beta;
    dwdx[5] = c;
}