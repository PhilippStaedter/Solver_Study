#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model22_karin1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = (1.0/10.0)*z;
    dwdx[1] = -mu*z;
    dwdx[2] = (1.0/10.0)*amici_y;
    dwdx[3] = lambdaMinus;
    dwdx[4] = -amici_y*mu;
}